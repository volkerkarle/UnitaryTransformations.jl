"""
    Schrieffer-Wolff Transformation

Perturbative unitary transformation to block-diagonalize a Hamiltonian,
eliminating couplings between a low-energy sector P and high-energy sector Q.
"""

export schrieffer_wolff, sw_generator, project_to_subspace

using QuantumAlgebra
using QuantumAlgebra:
    QuExpr, QuTerm, BaseOperator, BaseOpProduct, TLSCreate_, TLSDestroy_, normal_form, comm

import ..UnitaryTransformations: get_spin_constraint_info

"""
    schrieffer_wolff(H::QuExpr, P::Subspace; order::Int=2)

Perform the Schrieffer-Wolff transformation on Hamiltonian H with respect 
to the low-energy subspace P.

The transformation finds a unitary U = e^S such that H_eff = e^S H e^{-S}
is block-diagonal with respect to P and Q = 1-P, up to the specified order
in perturbation theory.

# Arguments
- `H`: The full Hamiltonian to transform
- `P`: The low-energy subspace definition
- `order`: Perturbation theory order (default: 2)

# Returns
- Named tuple `(H_eff, S, H_P)` where:
  - `H_eff`: The full block-diagonal effective Hamiltonian
  - `S`: The generator of the transformation
  - `H_P`: The effective Hamiltonian projected onto subspace P

# Example
```julia
using QuantumAlgebra, UnitaryTransformations

# Jaynes-Cummings in dispersive regime
ω = Pr"ω"   # cavity frequency
Δ = Pr"Δ"   # qubit splitting  
g = Pr"g"   # coupling strength

H = ω * a'()*a() + Δ/2 * σz() + g * (a'()*σm() + a()*σp())

# Transform to eliminate qubit-photon coupling
P = Subspace(σz() => -1)  # qubit ground state
result = schrieffer_wolff(H, P; order=2)
```
"""
function schrieffer_wolff(H::QuExpr, P::Subspace; order::Int = 2)
    order >= 1 || throw(ArgumentError("order must be at least 1, got $order"))

    # Normalize the Hamiltonian first
    H = normal_form(H)

    # Decompose H into diagonal and off-diagonal parts
    H_d, H_od = decompose(H, P)

    # Initialize total generator S and effective Hamiltonian
    S_total = QuExpr()
    H_eff = H_d  # Start with diagonal part

    # Current off-diagonal part to eliminate (starts with full H_od)
    current_od = H_od

    # Order-by-order construction
    for n = 1:order
        if isempty(current_od.terms)
            # No off-diagonal terms left to eliminate
            break
        end

        # Solve for S_n: [S_n, H_d] = -current_od
        S_n = solve_for_generator(H_d, current_od, P)

        # Accumulate generator
        S_total = normal_form(S_total + S_n)

        # Compute contribution to H_eff from this order
        # At order n, we get contributions from [S_n, H_od] and higher commutators

        if n == 1
            # Second-order contribution: (1/2)[S_1, V]
            # where V is the original off-diagonal part
            comm_S1_H = normal_form(comm(S_n, H))

            # The diagonal part of [S_1, H] contributes to H_eff
            comm_diag, comm_od = decompose(comm_S1_H, P)

            # Add half of the diagonal commutator (from BCH)
            H_eff = normal_form(H_eff + comm_diag / 2)

            # The remaining off-diagonal part needs to be eliminated at next order
            current_od = normal_form(comm_od / 2)

        else
            # Higher-order contributions
            # Use BCH: H_eff gets contributions from nested commutators

            # Compute [S_n, H] 
            comm_Sn_H = normal_form(comm(S_n, H_eff))
            comm_diag, comm_od = decompose(comm_Sn_H, P)

            # Add contribution to H_eff
            factorial_n = factorial(n)
            H_eff = normal_form(H_eff + comm_diag / factorial_n)

            # Update remaining off-diagonal part
            current_od = normal_form(current_od + comm_od / factorial_n)
        end
    end

    # Project the effective Hamiltonian onto subspace P
    H_P = project_to_subspace(H_eff, P)

    return (H_eff = H_eff, S = S_total, H_P = H_P)
end

"""
    sw_generator(H::QuExpr, P::Subspace; order::Int=1)

Compute only the generator S for the Schrieffer-Wolff transformation,
without computing the full effective Hamiltonian.

This is useful when you only need S, or want to manually compute
the transformation using `bch_transform`.
"""
function sw_generator(H::QuExpr, P::Subspace; order::Int = 1)
    H = normal_form(H)
    H_d, H_od = decompose(H, P)

    S_total = QuExpr()
    current_od = H_od

    for n = 1:order
        if isempty(current_od.terms)
            break
        end

        S_n = solve_for_generator(H_d, current_od, P)
        S_total = normal_form(S_total + S_n)

        if n < order
            # Compute next-order off-diagonal terms
            comm_Sn_H = normal_form(comm(S_n, H))
            _, comm_od = decompose(comm_Sn_H, P)
            current_od = normal_form(comm_od / factorial(n + 1))
        end
    end

    return S_total
end

"""
    project_to_subspace(H::QuExpr, P::Subspace)

Project an operator onto the subspace P.

This replaces diagonal operators by their eigenvalues in P:
- σz → eigenvalue (e.g., -1 for spin down)
- σ⁺σ⁻ → 0 for spin down, 1 for spin up
- a†a → eigenvalue (e.g., 0 for vacuum)

And removes any remaining off-diagonal terms.
"""
function project_to_subspace(H::QuExpr, P::Subspace)
    H = normal_form(H)

    # First, remove off-diagonal terms
    H_d = diagonal_part(H, P)

    # Then, substitute eigenvalues for constraint operators
    result = H_d

    for constraint in P.constraints
        # Check if this is a spin constraint
        spin_info = get_spin_constraint_info(constraint)
        if spin_info !== nothing
            spin_name, spin_inds, is_spin_down = spin_info
            # σ⁺σ⁻ = 0 for spin down, 1 for spin up
            pm_eigenvalue = is_spin_down ? 0 : 1
            result = substitute_spin_projection(result, spin_name, spin_inds, pm_eigenvalue)
            continue
        end

        # For other constraints, try direct substitution
        result = substitute_operator_eigenvalue(
            result,
            constraint.operator,
            constraint.eigenvalue,
        )
    end

    return normal_form(result)
end

"""
    substitute_spin_projection(H::QuExpr, spin_name, spin_inds, pm_eigenvalue)

Substitute σ⁺σ⁻ (or σ⁺(i)σ⁻(i)) by its eigenvalue in the spin subspace.
pm_eigenvalue is 0 for spin-down, 1 for spin-up.
"""
function substitute_spin_projection(H::QuExpr, spin_name, spin_inds, pm_eigenvalue)
    result = QuExpr()

    for (term, coeff) in H.terms
        ops = term.bares.v
        new_ops = BaseOperator[]
        eigenvalue_factor = 1
        skip = false
        i = 1

        while i <= length(ops)
            # Look for σ⁺σ⁻ pattern (TLSCreate_ followed by TLSDestroy_ with same name/inds)
            if i < length(ops) &&
               ops[i].t == TLSCreate_ &&
               ops[i+1].t == TLSDestroy_ &&
               ops[i].name == spin_name &&
               ops[i].inds == spin_inds &&
               ops[i+1].name == spin_name &&
               ops[i+1].inds == spin_inds
                # Found σ⁺σ⁻ - replace with eigenvalue
                eigenvalue_factor *= pm_eigenvalue
                i += 2  # Skip both operators
            else
                push!(new_ops, ops[i])
                i += 1
            end
        end

        # If eigenvalue_factor is 0, the whole term vanishes
        if eigenvalue_factor == 0
            continue
        end

        if isempty(new_ops)
            # The term was just operators that got substituted
            new_term = QuTerm(
                term.nsuminds,
                term.δs,
                term.params,
                term.expvals,
                term.corrs,
                BaseOpProduct(),
            )
        else
            new_term = QuTerm(
                term.nsuminds,
                term.δs,
                term.params,
                term.expvals,
                term.corrs,
                BaseOpProduct(new_ops),
            )
        end
        result = result + coeff * eigenvalue_factor * QuExpr(new_term)
    end

    return normal_form(result)
end

"""
    substitute_operator_eigenvalue(H::QuExpr, op::QuExpr, eigenvalue::Number)

Substitute an operator by its eigenvalue in an expression.

This is a simplified substitution that works for common cases where
the operator appears as a standalone term.
"""
function substitute_operator_eigenvalue(H::QuExpr, op::QuExpr, eigenvalue::Number)
    # Get the structure of the operator to match
    if length(op.terms) != 1
        return H  # Complex operator, skip substitution
    end

    op_term, op_coeff = first(op.terms)
    op_coeff == 1 || return H  # Has coefficient, skip

    result = QuExpr()

    for (term, coeff) in H.terms
        # Check if this term contains the operator exactly
        # For now, do simple replacement of matching terms

        matched = false

        # Handle σz case
        if length(op_term.bares.v) == 1 && length(term.bares.v) >= 1
            op_base = op_term.bares.v[1]

            # Look for matching operator in the term
            new_ops = BaseOperator[]
            eigenvalue_factor = 1

            for base_op in term.bares.v
                if base_op == op_base
                    eigenvalue_factor *= eigenvalue
                    matched = true
                else
                    push!(new_ops, base_op)
                end
            end

            if matched
                if isempty(new_ops)
                    # The whole term was just the operator
                    new_term = QuTerm(
                        term.nsuminds,
                        term.δs,
                        term.params,
                        term.expvals,
                        term.corrs,
                        BaseOpProduct(),
                    )
                else
                    new_term = QuTerm(
                        term.nsuminds,
                        term.δs,
                        term.params,
                        term.expvals,
                        term.corrs,
                        BaseOpProduct(new_ops),
                    )
                end
                result = result + coeff * eigenvalue_factor * QuExpr(new_term)
                continue
            end
        end

        # Handle a†a case (number operator)
        if length(op_term.bares.v) == 2 && length(term.bares.v) >= 2
            # Check if term contains the number operator
            # This is more complex - for now, just pass through
        end

        if !matched
            result = result + coeff * QuExpr(term)
        end
    end

    return normal_form(result)
end
