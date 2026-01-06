"""
    Schrieffer-Wolff Transformation

Perturbative unitary transformation to block-diagonalize a Hamiltonian,
eliminating couplings between a low-energy sector P and high-energy sector Q.
"""

export schrieffer_wolff, sw_generator, project_to_subspace

using QuantumAlgebra
using QuantumAlgebra:
    QuExpr, QuTerm, BaseOperator, BaseOpProduct, TLSCreate_, TLSDestroy_, normal_form, comm

import ..UnitaryTransformations:
    get_spin_constraint_info, simplify_coefficients, multi_nested_commutator, compositions

using Symbolics: Num, expand

"""
    _simplify_expand_only(expr::QuExpr)

Fast incremental simplification that only expands expressions.
This prevents expression tree explosion without expensive simplification.
"""
function _simplify_expand_only(expr::QuExpr)
    result_terms = Dict{QuTerm,Number}()
    for (term, coeff) in expr.terms
        if coeff isa Num
            # Just expand - this flattens nested sums/products efficiently
            simplified = expand(coeff)
            result_terms[term] = simplified
        else
            result_terms[term] = coeff
        end
    end
    return QuExpr(result_terms)
end

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
  Simplifying S can be slow at high orders. Set to `true` if you need simplified S.
- `simplify_mode`: Simplification mode for final output (default: `:fast`).
  - `:none` - No simplification (fastest)
  - `:fast` - Expansion only (default, very fast, recommended)
  - `:standard` - Basic algebraic simplification with expansion
  - `:fractions` - Simplify fractions with GCD (slow on complex expressions)
  - `:aggressive` - Full simplification (slowest)
- `diagonal_only`: If true, only compute H_eff (skip computing higher-order generators).
  This is much faster for high orders when you only need the effective Hamiltonian.
  Note: This uses a simplified algorithm that only uses S₁.
- `parallel`: If true, use multi-threading for orders > 3 (default: false).
  Requires Julia to be started with multiple threads (e.g., `julia -t 4`).
  For best performance, use 4-8 threads; more threads can cause lock contention.
  Parallelization is most beneficial for orders 4-6 with complex Hamiltonians.

# Returns
- Named tuple `(H_eff, S, H_P)` where:
  - `H_eff`: The full block-diagonal effective Hamiltonian
  - `S`: The generator of the transformation (S = S₁ + S₂ + ... where Sₙ is O(gⁿ))
        Note: If `diagonal_only=true`, S contains only S₁.
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

# For 4th order corrections (Kerr effect, etc.) - parallel computation
result4 = schrieffer_wolff(H, P; order=4, parallel=true)
```
"""
function schrieffer_wolff(
    H::QuExpr,
    P::Subspace;
    order::Int = 2,
    simplify_generator::Bool = false,
    simplify_mode::Symbol = :fast,
    diagonal_only::Bool = false,
    parallel::Bool = false,
)
    order >= 2 || throw(ArgumentError("order must be at least 2, got $order"))

    # Normalize the Hamiltonian first
    H = normal_form(H)

    # Decompose H = H₀ (diagonal) + V (off-diagonal)
    H0, V = decompose(H, P)

    # Use optimized diagonal-only algorithm if requested
    if diagonal_only
        return _schrieffer_wolff_diagonal_only(H, H0, V, P, order, simplify_mode)
    end

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

        order_n_terms = _collect_bch_terms_at_order(S, V_od, n, H0, V, P; parallel = parallel)
        # Simplify incrementally to prevent expression explosion
        # Use expand-only simplification which is fast but effective
        order_n_terms = _simplify_expand_only(order_n_terms)

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

    # Final simplification - only at the end, using specified mode
    H_eff = simplify_coefficients(H_eff; mode = simplify_mode)

    if simplify_generator
        S_total = simplify_coefficients(S_total; mode = simplify_mode)
    end

    # Project the effective Hamiltonian onto subspace P
    H_P = simplify_coefficients(project_to_subspace(H_eff, P); mode = simplify_mode)

    return (H_eff = H_eff, S = S_total, H_P = H_P)
end

"""
    _schrieffer_wolff_diagonal_only(H, H0, V, P, order, simplify_mode)

Optimized SW algorithm that only computes diagonal contributions to H_eff.

This is much faster for high orders because:
1. We only need S₁ (first-order generator)
2. We don't need to solve for higher-order generators
3. We only extract diagonal parts from nested commutators

The diagonal contributions at each order come from:
- Order 2: (1/2)[S₁, V] (diagonal part)
- Order 3: (1/6)[S₁, [S₁, V]] (diagonal part)  
- Order 4: (1/24)[S₁, [S₁, [S₁, V]]] (diagonal part)
- etc.

This uses the simplified BCH expansion with only S₁, which gives the leading
contributions at each order. Higher-order generators contribute at higher orders.
"""
function _schrieffer_wolff_diagonal_only(
    H::QuExpr,
    H0::QuExpr,
    V::QuExpr,
    P::Subspace,
    order::Int,
    simplify_mode::Symbol,
)
    # Compute only S₁
    S1 = solve_for_generator(H0, V, P)

    # Initialize effective Hamiltonian
    H_eff = H0

    # Compute nested commutators [S₁, [S₁, [..., [S₁, V]...]]] for increasing depth
    # and extract diagonal parts

    # Order 2: (1/2)[S₁, V]
    # Order 3: (1/6)[S₁, [S₁, V]] = (1/6)[S₁, C₁] where C₁ = [S₁, V]
    # Order k: (1/k!)[S₁, [S₁, [..., [S₁, V]...]]] (k-1 nested commutators with S₁)

    current_comm = V  # Start with V

    for k = 2:order
        # Compute [S₁, current_comm]
        next_comm = normal_form(comm(S1, current_comm))

        # BCH coefficient for this order
        bch_coeff = big(1) // factorial(big(k))

        # Extract diagonal part and add to H_eff
        diag_part, _ = decompose(bch_coeff * next_comm, P)
        H_eff = normal_form(H_eff + diag_part)

        # Update for next iteration
        current_comm = next_comm
    end

    # Final simplification
    H_eff = simplify_coefficients(H_eff; mode = simplify_mode)

    # Project to subspace
    H_P = simplify_coefficients(project_to_subspace(H_eff, P); mode = simplify_mode)

    return (H_eff = H_eff, S = S1, H_P = H_P)
end

"""
    WorkItem

A single unit of work for parallel BCH term computation.
Contains all information needed to compute one nested commutator term.
"""
struct WorkItem
    comp::Vector{Int}           # Composition (generator indices)
    base_type::Symbol           # :V or :H0
    bch_coeff::Rational{BigInt} # BCH coefficient (1/k!)
    inner_idx::Int              # For :H0 type, the index into V_od (0 for :V type)
end

"""
    _generate_work_items(S, V_od, n)

Pre-generate all work items for order n BCH term collection.
Returns a vector of WorkItem structs that can be processed in parallel.
"""
function _generate_work_items(S::Vector{QuExpr}, V_od::Vector{QuExpr}, n::Int)
    work_items = WorkItem[]
    max_depth = n

    for k = 1:max_depth
        bch_coeff = big(1) // factorial(big(k))

        # Part 1: Base operator is V (order 1)
        target_sum = n - 1
        if target_sum >= k
            for comp in compositions(target_sum, k; min_val = 1, max_val = n-1)
                if all(i -> isassigned(S, i), comp)
                    push!(work_items, WorkItem(comp, :V, bch_coeff, 0))
                end
            end
        end

        # Part 2: Base operator is H₀ (order 0)
        target_sum_H0 = n
        if target_sum_H0 >= k
            for comp in compositions(target_sum_H0, k; min_val = 1, max_val = n-1)
                if all(i -> isassigned(S, i), comp)
                    inner_idx = comp[end]
                    if isassigned(V_od, inner_idx) && !isempty(V_od[inner_idx].terms)
                        push!(work_items, WorkItem(comp, :H0, bch_coeff, inner_idx))
                    end
                end
            end
        end
    end

    return work_items
end

"""
    _compute_work_item(item, S, V_od, V)

Compute a single BCH term from a WorkItem.
This function is designed to be called from multiple threads.
"""
function _compute_work_item(
    item::WorkItem,
    S::Vector{QuExpr},
    V_od::Vector{QuExpr},
    V::QuExpr,
)
    if item.base_type == :V
        generators = QuExpr[S[i] for i in item.comp]
        term = multi_nested_commutator(generators, V)
        return item.bch_coeff * term
    else  # :H0
        inner_result = -V_od[item.inner_idx]
        if length(item.comp) == 1
            return item.bch_coeff * inner_result
        else
            outer_generators = QuExpr[S[i] for i in item.comp[1:(end-1)]]
            term = multi_nested_commutator(outer_generators, inner_result)
            return item.bch_coeff * term
        end
    end
end

"""
    _collect_bch_terms_at_order(S, V_od, n, H0, V, P; parallel=false)

Collect all BCH expansion terms at order n in the coupling.

The BCH expansion is: H_eff = H + [S,H] + (1/2)[S,[S,H]] + (1/6)[S,[S,[S,H]]] + ...

At order n, contributions come from nested commutators involving S₁,...,Sₙ₋₁
with appropriate BCH coefficients (1/k!).

This function dynamically computes all relevant terms by enumerating compositions
of the required order among the available generators.

# Arguments
- `parallel`: If true and n > 3, use multi-threading for term computation.

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
    P::Subspace;
    parallel::Bool = false,
)
    # Generate all work items
    work_items = _generate_work_items(S, V_od, n)

    # Use parallel execution for n > 3 when requested and multiple threads available
    # Note: For best performance with many threads (>8), consider limiting thread count
    # to avoid lock contention in QuantumAlgebra. Start Julia with: julia -t 4 or julia -t 8
    num_items = length(work_items)
    use_parallel = parallel && n > 3 && Threads.nthreads() > 1 && num_items > 1

    if use_parallel
        # Parallel execution: compute each term independently
        results = Vector{QuExpr}(undef, num_items)

        Threads.@threads for i in eachindex(work_items)
            results[i] = _compute_work_item(work_items[i], S, V_od, V)
        end

        # Serial reduction to avoid contention
        result = QuExpr()
        for r in results
            result = normal_form(result + r)
        end
        return result
    else
        # Sequential execution (original behavior)
        result = QuExpr()
        for item in work_items
            term = _compute_work_item(item, S, V_od, V)
            result = normal_form(result + term)
        end
        return result
    end
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
