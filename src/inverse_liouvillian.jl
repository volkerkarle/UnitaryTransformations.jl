"""
    Inverse Liouvillian: Solve [S, H_d] = -V_od for the generator S.

The key operation in Schrieffer-Wolff is finding S such that the commutator
with the diagonal Hamiltonian cancels the off-diagonal perturbation.

Two approaches are supported:
1. Eigenoperator method: For operators O where [H_d, O] = ε·O (TLS, bosons)
2. Matrix-element method: For general Lie algebras using ⟨i|L⁻¹[X]|j⟩ = Xᵢⱼ/(Eᵢ-Eⱼ)

Uses Symbolics.jl for proper handling of energy denominators like g²/(Δ-ω).
"""

export solve_for_generator,
    solve_for_generator_lie,
    solve_for_generator_eigenoperator,
    solve_for_generator_general,
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
    make_bare_quterm(bares::BaseOpProduct)

Create a QuTerm with only bare operators (no parameters, deltas, etc.).
This is a helper to avoid issues with internal QuantumAlgebra types.
"""
function make_bare_quterm(bares::BaseOpProduct)
    QuTerm(
        Int32(0),           # nsuminds
        QuantumAlgebra.δ[], # δs
        Param[],            # params
        ExpVal[],           # expvals
        Corr[],             # corrs
        bares,               # bares
    )
end

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

For an off-diagonal operator O, if [H_d, O] = ε·O, then the energy denominator is ε.
This is computed by evaluating [H_d, O] and extracting the coefficient.

Returns the energy denominator as a Symbolics Num expression.
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

    # Could not find matching term
    return nothing
end

"""
    solve_for_generator(H_d::QuExpr, V_od::QuExpr, P::Subspace; method::Symbol=:auto)

Solve [S, H_d] = -V_od for the generator S.

This is the core equation in Schrieffer-Wolff: we need to find S such that
the commutator with the diagonal part cancels the off-diagonal perturbation.

## Methods

Three methods are available, selected by the `method` keyword:

1. `:eigenoperator` - For operators O where [H_d, O] = ε·O (TLS, bosons).
   Solve S = O/ε directly for each term.

2. `:lie` - For SU(N) Lie algebras. Work in the Cartan-Weyl basis where 
   transition operators are eigenoperators: S_{ij} = V_{ij}/(E_i - E_j).

3. `:general` - For arbitrary operators. Construct the Liouvillian matrix
   L where [H_d, O_j] = Σ_i L_{ij} O_i and solve the linear system L·s = -v.

4. `:auto` (default) - Automatically select the best method:
   - Use `:lie` if Lie algebra generators are detected
   - Try `:eigenoperator` first, verify it works
   - Fall back to `:general` if eigenoperator method fails

Uses Symbolics.jl for proper symbolic division, allowing denominators like (Δ - ω).

# Arguments
- `H_d`: The diagonal (unperturbed) Hamiltonian
- `V_od`: The off-diagonal perturbation to be eliminated
- `P`: The subspace defining the block structure
- `method`: Solution method (:auto, :eigenoperator, :lie, :general)

# Returns
- `S`: The generator of the transformation
"""
function solve_for_generator(H_d::QuExpr, V_od::QuExpr, P::Subspace; method::Symbol = :auto)
    if method == :general
        return solve_for_generator_general(H_d, V_od, P)
    end

    if method == :eigenoperator
        return solve_for_generator_eigenoperator(H_d, V_od, P)
    end

    if method == :lie
        lie_info = detect_lie_algebra_system(V_od)
        if lie_info === nothing
            error("No Lie algebra operators detected in V_od, cannot use :lie method")
        end
        generators = get_generators_for_lie_system(lie_info)
        return solve_for_generator_lie(
            H_d,
            V_od,
            lie_info.N,
            generators;
            algebra_id = lie_info.algebra_id,
        )
    end

    # :auto method - try methods in order of preference

    # First, check if V_od contains Lie algebra operators
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

    # Try eigenoperator method first (fast path for TLS/bosons)
    S = solve_for_generator_eigenoperator(H_d, V_od, P)

    # Verify the solution: [S, H_d] should equal -V_od
    is_valid, residual = verify_generator(S, H_d, V_od)

    if is_valid
        return S
    end

    # Check if residual is small (some terms might just need simplification)
    residual_simplified = simplify_residual(residual)
    if isempty(residual_simplified.terms)
        return S
    end

    # Eigenoperator method failed - fall back to general method
    @debug "Eigenoperator method incomplete, using general method" residual=residual_simplified
    return solve_for_generator_general(H_d, V_od, P)
end

"""
    simplify_residual(residual::QuExpr)

Simplify the residual expression to check if it's effectively zero.
"""
function simplify_residual(residual::QuExpr)
    result_terms = Dict{QuTerm,Number}()
    for (term, coeff) in residual.terms
        simplified = Symbolics.simplify(coeff)
        if !iszero(simplified)
            result_terms[term] = simplified
        end
    end
    return QuExpr(result_terms)
end

"""
    solve_for_generator_eigenoperator(H_d::QuExpr, V_od::QuExpr, P::Subspace)

Solve [S, H_d] = -V_od using the eigenoperator method.

This method works for operators O where [H_d, O] = ε·O (eigenoperators of the adjoint action).
Examples include σ± for TLS, and bosonic operators a, a† when H_d is a number operator.

For each term in V_od, if [H_d, O] = ε·O, then the corresponding term in S is O/ε.
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

        # Compute energy denominator from [H_d, bare_O]
        denominator = compute_energy_denominator(H_d, term, P)

        if denominator === nothing
            bare_O = QuExpr(bare_term)
            @warn "Could not compute energy denominator for term: $bare_O"
            continue
        end

        # The coefficient for S is numerator / denominator
        # Don't simplify here to avoid expensive GCD computation
        # Simplification will happen at the end
        s_coeff = numerator / denominator

        # Add to result dictionary directly
        result_terms[bare_term] = s_coeff
    end

    return QuExpr(result_terms)
end

"""
    solve_for_generator_general(H_d::QuExpr, V_od::QuExpr, P::Subspace)

Solve [S, H_d] = -V_od for the generator S using the general method.

This method works for arbitrary operators where [H_d, O] may produce a linear 
combination of operators (not just a scalar multiple of O).

## Algorithm

1. Extract the operator basis {O_i} from V_od
2. Compute the Liouvillian matrix L where [H_d, O_j] = Σ_i L_{ij} O_i
3. Solve the linear system L·s = -v where v are coefficients of V_od in the basis

This reduces to the eigenoperator method when L is diagonal.

## When to use

- When operators mix under commutation with H_d (e.g., coupled modes)
- When the eigenoperator method fails to find matching terms in [H_d, O]
- As a fallback for complex operator structures

## Returns
- `S`: The generator of the transformation
"""
function solve_for_generator_general(H_d::QuExpr, V_od::QuExpr, P::Subspace)
    # Step 1: Extract operator basis from V_od
    # Each term in V_od defines a basis operator (the bare operator part)
    # We work with normal-ordered operators

    V_od = normal_form(V_od)

    if isempty(V_od.terms)
        return QuExpr()  # Nothing to solve
    end

    # Collect basis: map from bare operator structure to index
    basis_ops = Dict{BaseOpProduct,Int}()
    basis_list = BaseOpProduct[]  # Ordered list for indexing
    v_coeffs = Num[]  # Coefficients of V_od in the basis

    for (term, coeff) in V_od.terms
        bares = term.bares
        if !haskey(basis_ops, bares)
            push!(basis_list, bares)
            basis_ops[bares] = length(basis_list)
        end
        idx = basis_ops[bares]
        # Ensure v_coeffs is long enough
        while length(v_coeffs) < idx
            push!(v_coeffs, Num(0))
        end
        # Add the symbolic coefficient
        v_coeffs[idx] += symbolic_coefficient(term, coeff)
    end

    n = length(basis_list)

    # Step 2: Compute the Liouvillian matrix L
    # L[i,j] = coefficient of O_i in [H_d, O_j]

    L = Matrix{Num}(undef, n, n)
    fill!(L, Num(0))

    for j = 1:n
        # Create a pure operator from basis_list[j]
        bare_op = QuExpr(make_bare_quterm(basis_list[j]))

        # Compute [O_j, H_d] (note the order - we need [S, H_d] = -V_od)
        comm_result = normal_form(comm(bare_op, H_d))

        # Extract coefficients in the basis
        for (term, coeff) in comm_result.terms
            bares = term.bares
            if haskey(basis_ops, bares)
                i = basis_ops[bares]
                L[i, j] += symbolic_coefficient(term, coeff)
            else
                # The commutator produced an operator not in our basis
                # This can happen if the basis is incomplete
                # For now, we'll extend the basis
                push!(basis_list, bares)
                basis_ops[bares] = length(basis_list)
                i = length(basis_list)

                # Extend L and v_coeffs
                L_new = Matrix{Num}(undef, i, i)
                fill!(L_new, Num(0))
                L_new[1:n, 1:n] = L
                L = L_new
                L[i, j] = symbolic_coefficient(term, coeff)

                push!(v_coeffs, Num(0))
                n = i
            end
        end
    end

    # Step 3: Solve L·s = -v
    # For symbolic matrices, we use Cramer's rule or direct inversion
    # when the matrix is small, otherwise we need symbolic linear algebra

    s_coeffs = solve_symbolic_linear_system(L, -v_coeffs)

    if s_coeffs === nothing
        @warn "Could not solve linear system for generator (singular Liouvillian)"
        return QuExpr()
    end

    # Step 4: Reconstruct S from the solution
    result_terms = Dict{QuTerm,Number}()

    for i = 1:length(basis_list)
        if i > length(s_coeffs)
            break
        end
        s_i = s_coeffs[i]
        # Skip near-zero coefficients
        if iszero(Symbolics.simplify(s_i))
            continue
        end

        bare_term = make_bare_quterm(basis_list[i])
        result_terms[bare_term] = s_i
    end

    return QuExpr(result_terms)
end

"""
    solve_symbolic_linear_system(A::Matrix{Num}, b::Vector{Num})

Solve the linear system A·x = b symbolically.

For small systems (n ≤ 4), uses Cramer's rule.
For larger systems, attempts symbolic Gaussian elimination.

Returns `nothing` if the system is singular.
"""
function solve_symbolic_linear_system(A::Matrix{Num}, b::Vector{Num})
    n = size(A, 1)

    if n == 0
        return Num[]
    end

    if n == 1
        if iszero(Symbolics.simplify(A[1, 1]))
            return nothing
        end
        return [b[1] / A[1, 1]]
    end

    # For n ≤ 4, use Cramer's rule (explicit formulas)
    if n == 2
        return solve_2x2(A, b)
    elseif n == 3
        return solve_3x3(A, b)
    elseif n == 4
        return solve_4x4(A, b)
    end

    # For larger systems, use symbolic Gaussian elimination
    return solve_gaussian(A, b)
end

"""
    solve_2x2(A, b)

Solve 2×2 system using Cramer's rule.
"""
function solve_2x2(A::Matrix{Num}, b::Vector{Num})
    det_A = A[1, 1]*A[2, 2] - A[1, 2]*A[2, 1]
    det_A_simplified = Symbolics.simplify(det_A)

    if iszero(det_A_simplified)
        return nothing
    end

    x1 = (b[1]*A[2, 2] - b[2]*A[1, 2]) / det_A
    x2 = (A[1, 1]*b[2] - A[2, 1]*b[1]) / det_A

    return [x1, x2]
end

"""
    solve_3x3(A, b)

Solve 3×3 system using Cramer's rule.
"""
function solve_3x3(A::Matrix{Num}, b::Vector{Num})
    # Determinant using Sarrus' rule
    det_A =
        A[1, 1]*(A[2, 2]*A[3, 3] - A[2, 3]*A[3, 2]) -
        A[1, 2]*(A[2, 1]*A[3, 3] - A[2, 3]*A[3, 1]) +
        A[1, 3]*(A[2, 1]*A[3, 2] - A[2, 2]*A[3, 1])

    det_A_simplified = Symbolics.simplify(det_A)
    if iszero(det_A_simplified)
        return nothing
    end

    # Cramer's rule for each variable
    x = Vector{Num}(undef, 3)
    for i = 1:3
        A_i = copy(A)
        A_i[:, i] = b
        det_i =
            A_i[1, 1]*(A_i[2, 2]*A_i[3, 3] - A_i[2, 3]*A_i[3, 2]) -
            A_i[1, 2]*(A_i[2, 1]*A_i[3, 3] - A_i[2, 3]*A_i[3, 1]) +
            A_i[1, 3]*(A_i[2, 1]*A_i[3, 2] - A_i[2, 2]*A_i[3, 1])
        x[i] = det_i / det_A
    end

    return x
end

"""
    solve_4x4(A, b)

Solve 4×4 system using Cramer's rule.
"""
function solve_4x4(A::Matrix{Num}, b::Vector{Num})
    det_A = det_4x4(A)
    det_A_simplified = Symbolics.simplify(det_A)

    if iszero(det_A_simplified)
        return nothing
    end

    x = Vector{Num}(undef, 4)
    for i = 1:4
        A_i = copy(A)
        A_i[:, i] = b
        x[i] = det_4x4(A_i) / det_A
    end

    return x
end

"""
    det_4x4(A)

Compute determinant of 4×4 matrix using cofactor expansion.
"""
function det_4x4(A::Matrix{Num})
    # Expand along first row
    d = Num(0)
    for j = 1:4
        minor = [A[i, k] for i in 2:4, k in 1:4 if k != j]
        minor_det =
            minor[1, 1]*(minor[2, 2]*minor[3, 3] - minor[2, 3]*minor[3, 2]) -
            minor[1, 2]*(minor[2, 1]*minor[3, 3] - minor[2, 3]*minor[3, 1]) +
            minor[1, 3]*(minor[2, 1]*minor[3, 2] - minor[2, 2]*minor[3, 1])
        d += (-1)^(1+j) * A[1, j] * minor_det
    end
    return d
end

"""
    solve_gaussian(A, b)

Solve linear system using symbolic Gaussian elimination with partial pivoting.
"""
function solve_gaussian(A::Matrix{Num}, b::Vector{Num})
    n = size(A, 1)

    # Create augmented matrix
    aug = Matrix{Num}(undef, n, n+1)
    aug[:, 1:n] = A
    aug[:, n+1] = b

    # Forward elimination
    for k = 1:(n-1)
        # Find pivot (simplified symbolic check)
        pivot_row = k
        for i = (k+1):n
            # In symbolic case, we just use the first non-zero entry
            if iszero(Symbolics.simplify(aug[pivot_row, k])) &&
               !iszero(Symbolics.simplify(aug[i, k]))
                pivot_row = i
            end
        end

        # Swap rows if needed
        if pivot_row != k
            aug[k, :], aug[pivot_row, :] = aug[pivot_row, :], aug[k, :]
        end

        # Check for zero pivot
        if iszero(Symbolics.simplify(aug[k, k]))
            continue  # Try to continue (might still work)
        end

        # Eliminate below pivot
        for i = (k+1):n
            if !iszero(Symbolics.simplify(aug[i, k]))
                factor = aug[i, k] / aug[k, k]
                for j = k:(n+1)
                    aug[i, j] -= factor * aug[k, j]
                end
            end
        end
    end

    # Back substitution
    x = Vector{Num}(undef, n)
    for i = n:-1:1
        if iszero(Symbolics.simplify(aug[i, i]))
            return nothing  # Singular
        end
        x[i] = aug[i, n+1]
        for j = (i+1):n
            x[i] -= aug[i, j] * x[j]
        end
        x[i] /= aug[i, i]
    end

    return x
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

For SU(3), the Gell-Mann matrices λ₁...λ₆ are converted to raising/lowering operators:
- E₁₂ = (λ₁ + iλ₄)/2,  E₂₁ = (λ₁ - iλ₄)/2  (couple states 1↔2)
- E₁₃ = (λ₂ + iλ₅)/2,  E₃₁ = (λ₂ - iλ₅)/2  (couple states 1↔3)
- E₂₃ = (λ₃ + iλ₆)/2,  E₃₂ = (λ₃ - iλ₆)/2  (couple states 2↔3)

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

        # Extract off-diagonal elements
        for i = 1:N
            for j = 1:N
                if i != j && !iszero(gm[i, j])
                    # Matrix element ⟨i|λₖ|j⟩ contributes to E_{ij}
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

So for S = S[i,j] |i⟩⟨j| + S[j,i] |j⟩⟨i| (Hermitian):
    Real Gell-Mann coeff = S[i,j] + S[j,i]
    Imag Gell-Mann coeff = i(S[i,j] - S[j,i])
"""
function cartan_weyl_to_gellmann(
    transitions::Dict{Tuple{Int,Int},T},
    N::Int,
    generators::Tuple,
) where {T}
    # SU(3) specific: map transition pairs to Gell-Mann generator pairs
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
        S_up = get(transitions, pair, Num(0))           # S[i,j] where i<j (raising)
        S_down = get(transitions, (pair[2], pair[1]), Num(0))  # S[j,i] (lowering)

        # Get the Gell-Mann generator indices
        real_gen, imag_gen = transition_to_gen[pair]

        # |i⟩⟨j| = λ_real + i λ_imag  (for i < j)
        # |j⟩⟨i| = λ_real - i λ_imag
        #
        # S = S_up |i⟩⟨j| + S_down |j⟩⟨i|
        #   = S_up (λ_real + i λ_imag) + S_down (λ_real - i λ_imag)
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
1. Converting V_od to Cartan-Weyl basis (transition operators E_{ij})
2. Computing S_{ij} = V_{ij} / (Eᵢ - Eⱼ) for each transition
3. Converting S back to Gell-Mann basis

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
