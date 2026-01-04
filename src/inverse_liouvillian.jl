"""
    Inverse Liouvillian: Solve [S, H_d] = -V_od for the generator S.

The key operation in Schrieffer-Wolff is finding S such that the commutator
with the diagonal Hamiltonian cancels the off-diagonal perturbation.

Uses Symbolics.jl for proper handling of energy denominators like g²/(Δ-ω).
"""

export solve_for_generator,
    compute_energy_denominator, param_to_symbolic, symbolic_coefficient, clear_param_cache!

using QuantumAlgebra
using QuantumAlgebra:
    QuExpr, QuTerm, BaseOperator, BaseOpProduct, Param, normal_form, comm, QuOpName

using Symbolics
using Symbolics: Num, @variables

# Cache for parameter -> Symbolics variable mapping
const _param_cache = Dict{String,Num}()

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
"""
function symbolic_coefficient(term::QuTerm, coeff::Number)
    result = coeff isa Num ? coeff : Num(coeff)
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

    # Compute [H_d, bare_O]
    commutator = normal_form(comm(H_d, bare_O))

    if isempty(commutator.terms)
        return nothing  # O commutes with H_d (degenerate case)
    end

    # Find the term in commutator that matches the operator structure
    for (comm_term, comm_coeff) in commutator.terms
        if comm_term.bares == term.bares
            # Extract the full symbolic coefficient
            return symbolic_coefficient(comm_term, comm_coeff)
        end
    end

    # Could not find matching term
    return nothing
end

"""
    solve_for_generator(H_d::QuExpr, V_od::QuExpr, P::Subspace)

Solve [S, H_d] = -V_od for the generator S.

This is the core equation in Schrieffer-Wolff: we need to find S such that
the commutator with the diagonal part cancels the off-diagonal perturbation.

For each term in V_od, if [H_d, O] = ε·O, then the corresponding term in S is O/ε.
(Note: [S, H_d] = [O/ε, H_d] = (1/ε)[O, H_d] = -(1/ε)·ε·O = -O matches -V_od)

Uses Symbolics.jl for proper symbolic division, allowing denominators like (Δ - ω).

# Arguments
- `H_d`: The diagonal (unperturbed) Hamiltonian
- `V_od`: The off-diagonal perturbation to be eliminated
- `P`: The subspace defining the block structure

# Returns
- `S`: The generator of the transformation (anti-Hermitian)
"""
function solve_for_generator(H_d::QuExpr, V_od::QuExpr, P::Subspace)
    S = QuExpr()

    for (term, coeff) in V_od.terms
        # Create bare operator (strip params to compute commutator cleanly)
        bare_O = QuExpr(
            QuTerm(term.nsuminds, term.δs, Param[], term.expvals, term.corrs, term.bares),
        )

        # Get the full symbolic coefficient of this term in V_od
        numerator = symbolic_coefficient(term, coeff)

        # Compute energy denominator from [H_d, bare_O]
        denominator = compute_energy_denominator(H_d, term, P)

        if denominator === nothing
            @warn "Could not compute energy denominator for term: $bare_O"
            continue
        end

        # The coefficient for S is numerator / denominator
        # This is now proper symbolic division!
        s_coeff = numerator / denominator

        # Add to generator: s_coeff * bare_O
        S = S + s_coeff * bare_O
    end

    return normal_form(S)
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
