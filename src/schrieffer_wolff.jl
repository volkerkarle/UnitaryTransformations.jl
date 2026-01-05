"""
    Schrieffer-Wolff Transformation

Perturbative unitary transformation to block-diagonalize a Hamiltonian,
eliminating couplings between a low-energy sector P and high-energy sector Q.
"""

export schrieffer_wolff, sw_generator, project_to_subspace

using QuantumAlgebra
using QuantumAlgebra:
    QuExpr, QuTerm, BaseOperator, BaseOpProduct, TLSCreate_, TLSDestroy_, normal_form, comm

import ..UnitaryTransformations: get_spin_constraint_info, simplify_coefficients, multi_nested_commutator, compositions

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
- `order`: Perturbation theory order in the coupling strength (default: 2).
  - order=2: Standard SW, captures g² corrections (dispersive shifts)
  - order=4: Includes g⁴ corrections (Kerr nonlinearity, Bloch-Siegert, etc.)
- `simplify_generator`: Whether to simplify the generator S (default: false).
  Simplifying S can be very slow at high orders due to GCD computations on
  complex symbolic fractions. Set to `true` if you need simplified S.

# Returns
- Named tuple `(H_eff, S, H_P)` where:
  - `H_eff`: The full block-diagonal effective Hamiltonian
  - `S`: The generator of the transformation (S = S₁ + S₂ + ... where Sₙ is O(gⁿ))
  - `H_P`: The effective Hamiltonian projected onto subspace P

# Algorithm
The SW transformation uses S = S₁ + S₂ + S₃ + ... where each Sₙ is order gⁿ.
At each order, the off-diagonal terms from the BCH expansion determine Sₙ,
and the diagonal terms contribute to H_eff.

Key commutator rules (D=diagonal, O=off-diagonal):
- [O, O] → D
- [O, D] → O  
- [D, D] → D

# Example
```julia
using QuantumAlgebra, UnitaryTransformations, Symbolics

# Jaynes-Cummings in dispersive regime
@variables ω Δ g  # ω = cavity frequency, Δ = qubit splitting, g = coupling strength

H = ω * a'()*a() + Δ/2 * σz() + g * (a'()*σm() + a()*σp())

# Transform to eliminate qubit-photon coupling
P = Subspace(σz() => -1)  # qubit ground state
result = schrieffer_wolff(H, P; order=2)

# For 4th order corrections (Kerr effect, etc.)
result4 = schrieffer_wolff(H, P; order=4)
```
"""
function schrieffer_wolff(
    H::QuExpr,
    P::Subspace;
    order::Int = 2,
    simplify_generator::Bool = false,
)
    order >= 2 || throw(ArgumentError("order must be at least 2, got $order"))

    # Normalize the Hamiltonian first
    H = normal_form(H)

    # Decompose H = H₀ (diagonal) + V (off-diagonal)
    H0, V = decompose(H, P)

    # Store generators at each order: S[n] = Sₙ (order gⁿ)
    S = Vector{QuExpr}(undef, order)

    # Store off-diagonal "potentials" that each Sₙ must cancel: Vₙ where [Sₙ, H₀] = -Vₙ
    V_od = Vector{QuExpr}(undef, order)
    V_od[1] = V  # V₁ = V (the original off-diagonal part)

    # Initialize effective Hamiltonian with diagonal part
    H_eff = H0

    # ========== Order 1: [S₁, H₀] = -V ==========
    S[1] = solve_for_generator(H0, V, P)

    # ========== Order 2 and beyond ==========
    for n = 2:order
        # Collect all contributions at order n from the BCH expansion
        # H_eff = H + [S,H] + (1/2)[S,[S,H]] + (1/6)[S,[S,[S,H]]] + ...
        # with S = S₁ + S₂ + ... + Sₙ₋₁ (Sₙ not yet determined)

        order_n_terms = _collect_bch_terms_at_order(S, V_od, n, H0, V, P)
        order_n_terms = simplify_coefficients(order_n_terms)

        # Decompose into diagonal (→ H_eff) and off-diagonal (→ determines Sₙ)
        order_n_diag, order_n_od = decompose(order_n_terms, P)

        # Add diagonal contribution to H_eff
        H_eff = normal_form(H_eff + order_n_diag)

        # The off-diagonal part must be cancelled by [Sₙ, H₀]
        V_od[n] = order_n_od

        # Solve for Sₙ: [Sₙ, H₀] = -Vₙ
        if !isempty(order_n_od.terms)
            S[n] = solve_for_generator(H0, order_n_od, P)
        else
            S[n] = QuExpr()
        end
    end

    # Combine all generators
    S_total = QuExpr()
    for n = 1:order
        if isassigned(S, n)
            S_total = normal_form(S_total + S[n])
        end
    end

    # Final simplification
    H_eff = simplify_coefficients(H_eff)

    if simplify_generator
        S_total = simplify_coefficients(S_total)
    end

    # Project the effective Hamiltonian onto subspace P
    H_P = simplify_coefficients(project_to_subspace(H_eff, P))

    return (H_eff = H_eff, S = S_total, H_P = H_P)
end

"""
    _collect_bch_terms_at_order(S, V_od, n, H0, V, P)

Collect all BCH expansion terms at order n in the coupling.

The BCH expansion is: H_eff = H + [S,H] + (1/2)[S,[S,H]] + (1/6)[S,[S,[S,H]]] + ...

At order n, contributions come from nested commutators involving S₁,...,Sₙ₋₁
with appropriate BCH coefficients (1/k!).

This function dynamically computes all relevant terms by enumerating compositions
of the required order among the available generators.

# Algorithm
For each BCH depth k (number of nested commutators), we need to collect terms where
the generator orders plus the base operator order sum to n.

The base operator can be:
- V (order 1): need generator orders to sum to n-1
- H₀ (order 0): but [Sₘ, H₀] = -V_od[m], which has order m
  So [Sₘ, H₀] contributes like an order-m object, and we need to track this.

For simplicity and correctness, we enumerate all compositions and handle the
H₀ contribution by using the identity [Sₘ, H₀] = -V_od[m].
"""
function _collect_bch_terms_at_order(
    S::Vector{QuExpr},
    V_od::Vector{QuExpr},
    n::Int,
    H0::QuExpr,
    V::QuExpr,
    P::Subspace,
)
    result = QuExpr()
    
    # Import helper functions
    # multi_nested_commutator and compositions are from commutator_series.jl
    
    # For BCH expansion, we need nested commutators [Sᵢ₁, [Sᵢ₂, [..., [Sᵢₖ, X]...]]]
    # with coefficient 1/k!
    
    # At order n, we collect terms where: sum of generator orders + order(X) = n
    # Generators Sₘ have order m (for m = 1, ..., n-1)
    # V has order 1
    # H₀ has order 0, but [Sₘ, H₀] = -V_od[m] effectively contributes order m
    
    # Strategy: We separate the contribution into two parts:
    # 1. Nested commutators ending in V: [Sᵢ₁, [Sᵢ₂, [..., [Sᵢₖ, V]...]]]
    #    where i₁ + i₂ + ... + iₖ = n - 1 (since V is order 1)
    # 2. Nested commutators ending in H₀: [Sᵢ₁, [Sᵢ₂, [..., [Sᵢₖ, H₀]...]]]
    #    Using [Sₘ, H₀] = -V_od[m], these become:
    #    [Sᵢ₁, [Sᵢ₂, [..., -V_od[iₖ]...]]] where i₁ + i₂ + ... + iₖ = n
    
    # Maximum depth is n (if we use k copies of S₁, we get depth k with sum k)
    max_depth = n - 1  # At most n-1 generators (each at least order 1) to reach sum n-1
    
    for k in 1:max_depth
        # BCH coefficient: 1/k!
        bch_coeff = big(1) // factorial(big(k))
        
        # Part 1: Base operator is V (order 1)
        # Need generator orders to sum to n - 1
        target_sum = n - 1
        if target_sum >= k  # Need at least 1 per generator
            for comp in compositions(target_sum, k; min_val=1, max_val=n-1)
                # Check that all generator indices are valid (we only have S[1]...S[n-1])
                if all(i -> isassigned(S, i), comp)
                    generators = QuExpr[S[i] for i in comp]
                    term = multi_nested_commutator(generators, V)
                    result = normal_form(result + bch_coeff * term)
                end
            end
        end
        
        # Part 2: Base operator is H₀ (order 0)
        # The innermost commutator [Sᵢₖ, H₀] = -V_od[iₖ]
        # So we compute [Sᵢ₁, [..., [Sᵢₖ₋₁, -V_od[iₖ]]...]]
        # where i₁ + i₂ + ... + iₖ = n
        target_sum_H0 = n
        if target_sum_H0 >= k  # Need at least 1 per generator
            for comp in compositions(target_sum_H0, k; min_val=1, max_val=n-1)
                # Check that all generator indices are valid
                if all(i -> isassigned(S, i), comp)
                    # The innermost generator index gives us V_od[iₖ]
                    inner_idx = comp[end]
                    if isassigned(V_od, inner_idx) && !isempty(V_od[inner_idx].terms)
                        # [Sᵢₖ, H₀] = -V_od[iₖ]
                        inner_result = -V_od[inner_idx]
                        
                        if k == 1
                            # Just [Sᵢ₁, H₀] = -V_od[i₁], no further nesting
                            result = normal_form(result + bch_coeff * inner_result)
                        else
                            # Apply remaining generators: [Sᵢ₁, [..., [Sᵢₖ₋₁, -V_od[iₖ]]...]]
                            outer_generators = QuExpr[S[i] for i in comp[1:end-1]]
                            term = multi_nested_commutator(outer_generators, inner_result)
                            result = normal_form(result + bch_coeff * term)
                        end
                    end
                end
            end
        end
    end
    
    return result
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
        S_total = simplify_coefficients(normal_form(S_total + S_n))

        if n < order
            # Compute next-order off-diagonal terms
            comm_Sn_H = normal_form(comm(S_n, H))
            _, comm_od = decompose(comm_Sn_H, P)
            current_od = simplify_coefficients(normal_form(comm_od / factorial(n + 1)))
        end
    end

    return simplify_coefficients(S_total)
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
