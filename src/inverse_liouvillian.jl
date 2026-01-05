"""
    Inverse Liouvillian: Solve [S, H_d] = -V_od for the generator S.

The key operation in Schrieffer-Wolff is finding S such that the commutator
with the diagonal Hamiltonian cancels the off-diagonal perturbation.

The fundamental formula in the energy eigenbasis is:
    S_{ij} = V_{ij} / (E_i - E_j)

This is derived from [H_d, X] = A ⟹ (E_i - E_j)X_{ij} = A_{ij}.

Two approaches are supported depending on the operator type:
1. Eigenoperator method: For operators O where [H_d, O] = ε·O (TLS, bosons, N-level)
2. Matrix-element method: For SU(N) Lie algebras, convert to transition basis first

Uses Symbolics.jl for proper handling of energy denominators like g²/(Δ-ω).
"""

export solve_for_generator,
    solve_for_generator_lie,
    solve_for_generator_eigenoperator,
    compute_energy_denominator,
    compute_energy_eigenvalues,
    detect_lie_algebra_system,
    param_to_symbolic,
    symbolic_coefficient,
    clear_param_cache!

using QuantumAlgebra
using QuantumAlgebra:
    QuExpr,
    QuTerm,
    BaseOperator,
    BaseOpProduct,
    Param,
    ExpVal,
    Corr,
    normal_form,
    comm,
    QuOpName,
    LieAlgebraGen_,
    SU2_ALGEBRA_ID,
    SU3_ALGEBRA_ID

using Symbolics
using Symbolics: Num, @variables, simplify_fractions

# Cache for parameter -> Symbolics variable mapping
const _param_cache = Dict{String,Num}()

"""
    detect_lie_algebra_system(V_od::QuExpr)

Detect if V_od contains Lie algebra generators and return information about the algebra.

Returns `nothing` if no Lie algebra operators found, or a NamedTuple with:
- `N`: dimension of the algebra (2 for SU(2), 3 for SU(3))
- `algebra_id`: the UInt16 algebra identifier
- `name`: the generator name (e.g., :λ)
- `inds`: the indices tuple
"""
function detect_lie_algebra_system(V_od::QuExpr)
    for (term, _) in V_od.terms
        ops = term.bares.v
        for op in ops
            if op.t == LieAlgebraGen_
                # Found a Lie algebra generator
                algebra_id = op.algebra_id
                N =
                    algebra_id == SU2_ALGEBRA_ID ? 2 :
                    algebra_id == SU3_ALGEBRA_ID ? 3 :
                    error("Unknown algebra_id: $algebra_id")
                return (N = N, algebra_id = algebra_id, name = op.name, inds = op.inds)
            end
        end
    end
    return nothing
end

"""
    get_generators_for_lie_system(lie_info::NamedTuple)

Get the SU(N) generators tuple for the detected Lie algebra system.
"""
function get_generators_for_lie_system(lie_info::NamedTuple)
    # Convert QuOpName to Symbol for su_generators
    name_sym = Symbol(lie_info.name)
    return su_generators(lie_info.N, name_sym)
end

"""
    param_to_symbolic(p::Param)

Convert a QuantumAlgebra Param to a Symbolics variable.
Caches variables to ensure the same param always maps to the same variable.
"""
function param_to_symbolic(p::Param)
    name_str = string(p.name)
    if haskey(_param_cache, name_str)
        return _param_cache[name_str]
    else
        var = Symbolics.variable(Symbol(name_str))
        _param_cache[name_str] = var
        return var
    end
end

"""
    symbolic_coefficient(term::QuTerm, coeff::Number)

Convert a QuTerm's coefficient (Number + Params) to a single Symbolics expression.
Handles both real and complex coefficients.
"""
function symbolic_coefficient(term::QuTerm, coeff::Number)
    # Handle different coefficient types
    if coeff isa Num
        result = coeff
    elseif coeff isa Complex
        # Handle complex coefficients - keep as Complex{Num} or Complex{Real}
        # For Complex{Num}, we need to preserve the structure
        if real(coeff) isa Num || imag(coeff) isa Num
            result = coeff  # Already has symbolic components
        else
            # Pure numeric complex - convert to Num
            result = Num(real(coeff)) + im * Num(imag(coeff))
        end
    else
        result = Num(coeff)
    end

    for p in term.params
        result = result * param_to_symbolic(p)
    end
    return result
end

"""
    compute_energy_denominator(H_d::QuExpr, term::QuTerm, P::Subspace)

Compute the energy denominator for a given off-diagonal term.

For an off-diagonal operator O, computes ε where [H_d, O] = ε·O.
This is the energy difference E_i - E_j for the transition that O represents.

Returns the energy denominator as a Symbolics Num expression, or `nothing` if
the operator commutes with H_d (degenerate case).
"""
function compute_energy_denominator(H_d::QuExpr, term::QuTerm, P::Subspace)
    # Create bare operator (strip params)
    bare_O = QuExpr(
        QuTerm(term.nsuminds, term.δs, Param[], term.expvals, term.corrs, term.bares),
    )

    # Compute [H_d, bare_O] - this will normal-order the result
    commutator = normal_form(comm(H_d, bare_O))

    if isempty(commutator.terms)
        return nothing  # O commutes with H_d (degenerate case)
    end

    # Normal-order the original term for comparison
    # (the commutator result is already normal-ordered)
    bare_O_normal = normal_form(bare_O)

    # Get the bares from the normal-ordered original
    if isempty(bare_O_normal.terms)
        return nothing
    end
    expected_bares = first(bare_O_normal.terms)[1].bares

    # Find the term in commutator that matches the operator structure
    for (comm_term, comm_coeff) in commutator.terms
        if comm_term.bares == expected_bares
            # Extract the full symbolic coefficient
            return symbolic_coefficient(comm_term, comm_coeff)
        end
    end

    # Could not find matching term - operator is not an eigenoperator
    return nothing
end

"""
    solve_for_generator(H_d::QuExpr, V_od::QuExpr, P::Subspace)

Solve [S, H_d] = -V_od for the generator S.

This is the core equation in Schrieffer-Wolff: we need to find S such that
the commutator with the diagonal part cancels the off-diagonal perturbation.

The fundamental solution in the energy eigenbasis is:
    S_{ij} = V_{ij} / (E_i - E_j)

Two methods are automatically selected based on operator types:

1. **Eigenoperator method** (TLS, bosons, N-level transitions): For operators O 
   where [H_d, O] = ε·O, the solution is S = V/ε directly.

2. **Matrix-element method** (SU(N) Lie algebras): Convert to the Cartan-Weyl 
   (transition) basis where operators ARE eigenoperators, apply S_{ij} = V_{ij}/(E_i - E_j),
   then convert back.

Uses Symbolics.jl for proper symbolic division, allowing denominators like (Δ - ω).

# Arguments
- `H_d`: The diagonal (unperturbed) Hamiltonian
- `V_od`: The off-diagonal perturbation to be eliminated
- `P`: The subspace defining the block structure

# Returns
- `S`: The generator of the transformation
"""
function solve_for_generator(H_d::QuExpr, V_od::QuExpr, P::Subspace)
    # Check if V_od contains Lie algebra operators
    lie_info = detect_lie_algebra_system(V_od)

    if lie_info !== nothing
        # Use matrix-element method for Lie algebras
        generators = get_generators_for_lie_system(lie_info)
        return solve_for_generator_lie(
            H_d,
            V_od,
            lie_info.N,
            generators;
            algebra_id = lie_info.algebra_id,
        )
    end

    # Use eigenoperator method for TLS/bosons/N-level transitions
    return solve_for_generator_eigenoperator(H_d, V_od, P)
end

"""
    solve_for_generator_eigenoperator(H_d::QuExpr, V_od::QuExpr, P::Subspace)

Solve [S, H_d] = -V_od using the eigenoperator method.

This method works for operators O where [H_d, O] = ε·O (eigenoperators of the 
adjoint action). The solution is simply S = V/ε.

Examples of eigenoperators:
- σ± for TLS with H_d ∝ σz: [σz, σ±] = ±2σ±
- a, a† for bosons with H_d ∝ a†a: [a†a, a] = -a, [a†a, a†] = a†
- |i⟩⟨j| for N-level with diagonal H_d: [H_d, |i⟩⟨j|] = (E_i - E_j)|i⟩⟨j|

For each term in V_od, the corresponding term in S is V_term / ε.
"""
function solve_for_generator_eigenoperator(H_d::QuExpr, V_od::QuExpr, P::Subspace)
    # Build result directly to avoid expensive iszero checks in + operator
    result_terms = Dict{QuTerm,Number}()

    for (term, coeff) in V_od.terms
        # Create bare operator (strip params to compute commutator cleanly)
        bare_term =
            QuTerm(term.nsuminds, term.δs, Param[], term.expvals, term.corrs, term.bares)

        # Get the full symbolic coefficient of this term in V_od
        numerator = symbolic_coefficient(term, coeff)

        # Compute energy denominator from [H_d, bare_O] = ε·bare_O
        denominator = compute_energy_denominator(H_d, term, P)

        if denominator === nothing
            bare_O = QuExpr(bare_term)
            @warn "Could not compute energy denominator for term: $bare_O"
            continue
        end

        # The coefficient for S is V / ε = numerator / denominator
        s_coeff = numerator / denominator

        # Add to result dictionary directly
        result_terms[bare_term] = s_coeff
    end

    return QuExpr(result_terms)
end

"""
    verify_generator(S::QuExpr, H_d::QuExpr, V_od::QuExpr)

Verify that the generator S satisfies [S, H_d] ≈ -V_od.

Returns true if the equation is satisfied (up to normal ordering), false otherwise.
Also returns the residual [S, H_d] + V_od which should be zero.
"""
function verify_generator(S::QuExpr, H_d::QuExpr, V_od::QuExpr)
    lhs = normal_form(comm(S, H_d))
    residual = normal_form(lhs + V_od)
    is_zero = isempty(residual.terms)
    return is_zero, residual
end

"""
    clear_param_cache!()

Clear the parameter to Symbolics variable cache.
Useful when starting a new calculation with fresh variables.
"""
function clear_param_cache!()
    empty!(_param_cache)
end

# =============================================================================
# Matrix-Element Method for Lie Algebras (Cartan-Weyl Basis)
# =============================================================================

"""
    compute_energy_eigenvalues(H_d::QuExpr, N::Int, algebra_id::UInt16)

Compute energy eigenvalues for an N-level system from a diagonal Hamiltonian
expressed in terms of Lie algebra generators.

For SU(N), the diagonal generators (Cartan subalgebra) have known eigenvalues.
This function extracts the coefficients from H_d and computes Eᵢ for each state.

Returns a vector of symbolic expressions [E₁, E₂, ..., Eₙ].
"""
function compute_energy_eigenvalues(H_d::QuExpr, N::Int, algebra_id::UInt16)
    eigenvalues = [Num(0) for _ = 1:N]

    for (term, coeff) in H_d.terms
        ops = term.bares.v

        # Skip identity terms (no operators)
        if isempty(ops)
            sym_coeff = symbolic_coefficient(term, coeff)
            for i = 1:N
                eigenvalues[i] += sym_coeff
            end
            continue
        end

        # Single Lie algebra generator term
        if length(ops) == 1 && ops[1].t == LieAlgebraGen_
            op = ops[1]
            if op.algebra_id == algebra_id
                gen_idx = Int(op.gen_idx)
                sym_coeff = symbolic_coefficient(term, coeff)

                # Get the diagonal elements of this generator
                gm = QuantumAlgebra.gellmann_matrix(N, gen_idx)
                for i = 1:N
                    eigenvalues[i] += sym_coeff * real(gm[i, i])
                end
            end
        end
    end

    return eigenvalues
end

"""
    gellmann_to_cartan_weyl(V_od::QuExpr, N::Int, algebra_id::UInt16)

Convert an off-diagonal operator from Gell-Mann basis to Cartan-Weyl basis.

The Cartan-Weyl basis uses transition operators |i⟩⟨j| which ARE eigenoperators
of the diagonal Hamiltonian: [H_d, |i⟩⟨j|] = (E_i - E_j)|i⟩⟨j|.

For SU(3), the Gell-Mann matrices λ₁...λ₆ are converted to:
- E₁₂ = |1⟩⟨2|, E₂₁ = |2⟩⟨1| (from λ₁, λ₄)
- E₁₃ = |1⟩⟨3|, E₃₁ = |3⟩⟨1| (from λ₂, λ₅)
- E₂₃ = |2⟩⟨3|, E₃₂ = |3⟩⟨2| (from λ₃, λ₆)

Returns a Dict mapping (i,j) => coefficient for transition i→j.
"""
function gellmann_to_cartan_weyl(V_od::QuExpr, N::Int, algebra_id::UInt16)
    # Map (i,j) -> symbolic coefficient
    transitions = Dict{Tuple{Int,Int},Any}()

    for (term, coeff) in V_od.terms
        ops = term.bares.v

        if length(ops) != 1 || ops[1].t != LieAlgebraGen_
            @warn "Non-Lie algebra term in V_od, skipping: $term"
            continue
        end

        op = ops[1]
        if op.algebra_id != algebra_id
            continue
        end

        gen_idx = Int(op.gen_idx)
        sym_coeff = symbolic_coefficient(term, coeff)

        # Get the matrix representation
        gm = QuantumAlgebra.gellmann_matrix(N, gen_idx)

        # Extract off-diagonal elements: V_{ij} = coeff * ⟨i|λₖ|j⟩
        for i = 1:N
            for j = 1:N
                if i != j && !iszero(gm[i, j])
                    key = (i, j)
                    matrix_elem = gm[i, j]
                    contribution = sym_coeff * matrix_elem

                    if haskey(transitions, key)
                        transitions[key] += contribution
                    else
                        transitions[key] = contribution
                    end
                end
            end
        end
    end

    return transitions
end

"""
    cartan_weyl_to_gellmann(transitions::Dict, N::Int, generators::Tuple)

Convert from Cartan-Weyl (transition) basis back to Gell-Mann basis.

For an operator S with matrix elements S[i,j] (transition from j to i),
we use the relation between matrix units and Gell-Mann matrices.

For SU(3), the off-diagonal matrix units are:
    |1⟩⟨2| = λ₁ + iλ₄,  |2⟩⟨1| = λ₁ - iλ₄
    |1⟩⟨3| = λ₂ + iλ₅,  |3⟩⟨1| = λ₂ - iλ₅
    |2⟩⟨3| = λ₃ + iλ₆,  |3⟩⟨2| = λ₃ - iλ₆
"""
function cartan_weyl_to_gellmann(
    transitions::Dict{Tuple{Int,Int},T},
    N::Int,
    generators::Tuple,
) where {T}
    # SU(N) specific: map transition pairs to Gell-Mann generator pairs
    # (i,j) with i<j maps to (real_gen, imag_gen)
    if N == 3
        transition_to_gen = Dict(
            (1, 2) => (1, 4),  # |1⟩⟨2| uses λ₁ + iλ₄
            (1, 3) => (2, 5),  # |1⟩⟨3| uses λ₂ + iλ₅
            (2, 3) => (3, 6),  # |2⟩⟨3| uses λ₃ + iλ₆
        )
    elseif N == 2
        transition_to_gen = Dict(
            (1, 2) => (1, 2),  # |1⟩⟨2| uses σ₁ + iσ₂
        )
    else
        error("cartan_weyl_to_gellmann not implemented for N=$N")
    end

    result = zero(QuExpr)
    processed = Set{Tuple{Int,Int}}()

    for ((i, j), S_ij) in transitions
        # Process each pair (i,j) and (j,i) together
        pair = i < j ? (i, j) : (j, i)
        if pair in processed
            continue
        end
        push!(processed, pair)

        # Get both transitions
        S_up = get(transitions, pair, Num(0))           # S[i,j] where i<j
        S_down = get(transitions, (pair[2], pair[1]), Num(0))  # S[j,i]

        # Get the Gell-Mann generator indices
        real_gen, imag_gen = transition_to_gen[pair]

        # |i⟩⟨j| = λ_real + i λ_imag  (for i < j)
        # |j⟩⟨i| = λ_real - i λ_imag
        #
        # S = S_up |i⟩⟨j| + S_down |j⟩⟨i|
        #   = (S_up + S_down) λ_real + i(S_up - S_down) λ_imag

        coeff_real = S_up + S_down
        coeff_imag = 1im * (S_up - S_down)

        result = result + coeff_real * generators[real_gen]
        result = result + coeff_imag * generators[imag_gen]
    end

    return normal_form(result)
end

"""
    solve_for_generator_lie(H_d::QuExpr, V_od::QuExpr, N::Int, generators::Tuple; 
                            algebra_id=SU3_ALGEBRA_ID)

Solve [S, H_d] = -V_od for the generator S using the matrix-element method.

This approach works for any SU(N) Lie algebra by:
1. Computing energy eigenvalues E_i from H_d
2. Converting V_od to Cartan-Weyl basis (transition operators |i⟩⟨j|)
3. Applying the fundamental formula: S_{ij} = V_{ij} / (E_i - E_j)
4. Converting S back to Gell-Mann basis

The key insight is that transition operators |i⟩⟨j| ARE eigenoperators:
    [H_d, |i⟩⟨j|] = (E_i - E_j)|i⟩⟨j|

# Arguments
- `H_d`: The diagonal Hamiltonian (in diagonal generators only)
- `V_od`: The off-diagonal perturbation to be eliminated
- `N`: Dimension of the representation (2 for SU(2), 3 for SU(3))
- `generators`: Tuple of SU(N) generators (from su_generators)
- `algebra_id`: The algebra identifier (default: SU3_ALGEBRA_ID)

# Returns
- `S`: The generator of the transformation
"""
function solve_for_generator_lie(
    H_d::QuExpr,
    V_od::QuExpr,
    N::Int,
    generators::Tuple;
    algebra_id::UInt16 = SU3_ALGEBRA_ID,
)
    # Step 1: Compute energy eigenvalues
    E = compute_energy_eigenvalues(H_d, N, algebra_id)

    # Step 2: Convert V_od to transition basis
    V_transitions = gellmann_to_cartan_weyl(V_od, N, algebra_id)

    # Step 3: Apply inverse Liouvillian in transition basis
    # S_{ij} = V_{ij} / (E_i - E_j)
    S_transitions = Dict{Tuple{Int,Int},Any}()

    for ((i, j), V_ij) in V_transitions
        denominator = E[i] - E[j]
        S_transitions[(i, j)] = V_ij / denominator
    end

    # Step 4: Convert S back to Gell-Mann basis
    S = cartan_weyl_to_gellmann(S_transitions, N, generators)

    return S
end
