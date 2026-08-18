"""
    Schrieffer-Wolff Transformation

Perturbative unitary transformation to block-diagonalize a Hamiltonian,
eliminating couplings between a low-energy sector P and high-energy sector Q.
"""

export schrieffer_wolff, sw_generator, project_to_subspace

using QuantumAlgebra
using QuantumAlgebra:
    QuExpr,
    QuTerm,
    BaseOperator,
    BaseOpProduct,
    TLSCreate_,
    TLSDestroy_,
    Transition_,
    normal_form,
    comm

import ..UnitaryTransformations:
    get_spin_constraint_info,
    get_transition_constraint_info,
    get_transition_indices,
    is_number_constraint,
    simplify_coefficients,
    multi_nested_commutator,
    compositions,
    is_fastmode_flat

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
    _transition_constraint_infos(P::Subspace)

Collect transition constraint info for include_QQ handling.
"""
function _transition_constraint_infos(P::Subspace)
    infos = NamedTuple[]
    for constraint in P.constraints
        info = get_transition_constraint_info(constraint)
        if info !== nothing
            push!(infos, info)
        end
    end
    return infos
end

"""
    _term_has_QQ_transition(term::QuTerm, transition_infos)

Return true if a term contains an off-diagonal transition operator that does
not involve the constrained state (eigenvalue == 1), i.e., a Q↔Q coupling.
"""
function _term_has_QQ_transition(term::QuTerm, transition_infos)
    isempty(transition_infos) && return false
    for op in term.bares.v
        op.t == Transition_ || continue
        indices = get_transition_indices(op)
        indices === nothing && continue
        i, j = indices
        i == j && continue
        for info in transition_infos
            if op.name == info.name && op.inds == info.inds
                if info.eigenvalue == 1 && i != info.state && j != info.state
                    return true
                end
            end
        end
    end
    return false
end

"""
    _extract_QQ_terms(H_d::QuExpr, P::Subspace)

Extract block-diagonal terms that couple within Q for transition constraints.
"""
function _extract_QQ_terms(H_d::QuExpr, P::Subspace)
    transition_infos = _transition_constraint_infos(P)
    isempty(transition_infos) && return QuExpr()
    result_terms = Dict{QuTerm,Number}()
    for (term, coeff) in H_d.terms
        if _term_has_QQ_transition(term, transition_infos)
            result_terms[term] = coeff
        end
    end
    return QuExpr(result_terms)
end

"""
    schrieffer_wolff(H::QuExpr, P::Subspace; order::Int=2, include_QQ::Bool=true, ...)

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
   Note: This uses a simplified algorithm that only uses S₁. The `include_QQ`
   parameter is ignored when `diagonal_only=true`; Q↔Q virtual paths are not included.
- `include_QQ`: If true (default), treat block-diagonal couplings within Q (e.g.,
  off-diagonal transition operators that do not involve the constrained state)
  as perturbations. This captures virtual paths like P→Q→Q→P at higher order.
  Set to false to use the legacy behavior (expand only in P↔Q couplings).
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
    include_QQ::Bool = true,
    parallel::Bool = false,
    fastmode_flat::Bool = false,
)
    order >= 2 || throw(ArgumentError("order must be at least 2, got $order"))
    if fastmode_flat && parallel
        @warn "fastmode_flat disables threaded BCH collection for Symbolics stability"
        parallel = false
    end

    prev_fast = is_fastmode_flat()
    set_fastmode_flat!(fastmode_flat)
    try
        # Normalize the Hamiltonian first
        H = normal_form(H)

        # Decompose H = H_d (block-diagonal) + H_od (P↔Q couplings)
        H_d, H_od = decompose(H, P)

        # Use optimized diagonal-only algorithm if requested
        if diagonal_only
            if include_QQ
                @warn "diagonal_only=true ignores include_QQ; Q↔Q virtual paths are not included"
            end
            return _schrieffer_wolff_diagonal_only(H, H_d, H_od, P, order, simplify_mode)
        end

        # Optionally treat Q↔Q couplings as perturbations
        H0 = H_d
        V_total = H_od
        if include_QQ
            V_QQ = _extract_QQ_terms(H_d, P)
            if !isempty(V_QQ.terms)
                H0 = normal_form(H_d - V_QQ)
                V_total = normal_form(H_od + V_QQ)
            end
        end

        # Split perturbation into diagonal and off-diagonal parts
        V_diag, V_od1 = decompose(V_total, P)

        # Store generators at each order: S[n] = Sₙ (order gⁿ)
        S = Vector{QuExpr}(undef, order)

        # Store off-diagonal "potentials" that each Sₙ must cancel: Vₙ where [Sₙ, H₀] = -Vₙ
        V_od = Vector{QuExpr}(undef, order)
        V_od[1] = V_od1  # V₁ = off-diagonal part of perturbation

        # Initialize effective Hamiltonian with diagonal part — accumulate as dict to avoid
        # expensive Symbolics polynomial simplification on every normal_form call
        H_eff_dict = Dict{QuTerm,Number}()
        for (term, coeff) in normal_form(H0 + V_diag).terms
            H_eff_dict[term] = coeff
        end

        # ========== Order 1: [S₁, H₀] = -V ==========
        S[1] = solve_for_generator(H0, V_od1, P)

        # ========== Order 2 and beyond ==========
        for n = 2:order
            # Collect all contributions at order n from the BCH expansion
            order_n_terms =
                _collect_bch_terms_at_order_sequential(S, V_od, n, H0, V_total, P; parallel = parallel)
            # Clear cache between orders to avoid unbounded growth
            _clear_bch_cache!()

            # NOTE: skip _simplify_expand_only to avoid Symbolics polynomial explosion.
            # Decompose works correctly on raw expressions.

            # Decompose into diagonal (→ H_eff) and off-diagonal (→ determines Sₙ)
            order_n_diag, order_n_od = decompose(order_n_terms, P)

            # Add diagonal contribution to H_eff dict (no normal_form)
            for (term, coeff) in order_n_diag.terms
                H_eff_dict[term] = get(H_eff_dict, term, 0) + coeff
            end

            # The off-diagonal part must be cancelled by [Sₙ, H₀]
            V_od[n] = order_n_od

            # Solve for Sₙ: [Sₙ, H₀] = -Vₙ
            if !isempty(order_n_od.terms)
                S[n] = solve_for_generator(H0, order_n_od, P)
            else
                S[n] = QuExpr()
            end
        end

        # Convert accumulated dict back to QuExpr (no normal_form to avoid Symbolics explosion)
        H_eff = QuExpr(H_eff_dict)

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
    finally
        set_fastmode_flat!(prev_fast)
    end
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
            for comp in compositions(target_sum, k; min_val = 1, max_val = n - 1)
                if all(i -> isassigned(S, i), comp)
                    push!(work_items, WorkItem(comp, :V, bch_coeff, 0))
                end
            end
        end

        # Part 2: Base operator is H₀ (order 0)
        target_sum_H0 = n
        if target_sum_H0 >= k
            for comp in compositions(target_sum_H0, k; min_val = 1, max_val = n - 1)
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

# ======================================================================
# Sequential BCH term collection with memoization (replaces compositions-based approach)
# Instead of enumerating all integer compositions (which explodes at order 4+),
# we compute nested commutators recursively and cache intermediate results.
# ======================================================================

# Memoization cache: (generator_indices_tuple, base_string) => QuExpr
const _bch_cache = Dict{Any,QuExpr}()
# Lock for thread safety
const _bch_cache_lock = ReentrantLock()

function _collapse_transition_ops(term::QuTerm)
    ops = term.bares.v
    isempty(ops) && return term

    # Transition operators commute with bosonic/TLS operators; group transitions
    # to expose products like Lᵢⱼ Lⱼₖ that can be reduced.
    nontrans = BaseOperator[]
    trans = BaseOperator[]
    for op in ops
        if op.t == Transition_
            push!(trans, op)
        else
            push!(nontrans, op)
        end
    end

    # Fast zero-prune: if two adjacent transition factors cannot compose
    # (L[a,b]L[c,d] with b!=c), the whole term is zero.
    for i = 1:(length(trans)-1)
        op1 = trans[i]
        op2 = trans[i + 1]
        idx1 = get_transition_indices(op1)
        idx2 = get_transition_indices(op2)
        if idx1 !== nothing && idx2 !== nothing && op1.name == op2.name && op1.inds == op2.inds
            _, b = idx1
            c, _ = idx2
            if b != c
                return nothing
            end
        end
    end

    new_ops = vcat(nontrans, trans)
    new_bares = isempty(new_ops) ? BaseOpProduct() : BaseOpProduct(new_ops)
    return QuTerm(term.nsuminds, term.δs, term.params, term.expvals, term.corrs, new_bares)
end

function _collapse_transition_terms(expr::QuExpr)
    result_terms = Dict{QuTerm,Number}()
    for (term, coeff) in expr.terms
        new_term = _collapse_transition_ops(term)
        new_term === nothing && continue
        result_terms[new_term] = get(result_terms, new_term, 0) + coeff
    end
    return QuExpr(result_terms)
end

function _clear_bch_cache!()
    lock(_bch_cache_lock) do
        empty!(_bch_cache)
    end
end

"""
    _nested_comm_memoized(gen_is::NTuple{K,Int}, base_symbol::String, S, V, V_od) where K

Compute the nested commutator [S_{i₁}, [S_{i₂}, ..., [S_{iₖ}, base]...]] where
base is V (if base_symbol = "V") or -V_od[inner_idx] (if base_symbol starts with "H0").

Results are cached in _bch_cache for reuse across orders.
"""
function _nested_comm_memoized(
    gen_is::NTuple{K,Int},
    base_symbol::String,
    S::Vector{QuExpr},
    V::QuExpr,
    V_od::Vector{QuExpr},
) where {K}
    # Check cache
    cache_key = (gen_is, base_symbol)
    cached = lock(_bch_cache_lock) do
        get(_bch_cache, cache_key, nothing)
    end
    if cached !== nothing
        return cached
    end

    if base_symbol == "V"
        # V-based: all generator indices go into the nested commutator against V
        result = multi_nested_commutator(QuExpr[S[i] for i in gen_is], V)
    else  # "H0:3" means innermost = -V_od[3], outer generators are gen_is[1:end-1]
        inner_idx = parse(Int, split(base_symbol, ":")[2])
        if K == 1
            # [S_i, H0] = -V_od[i] — just return the base directly (no outer generators)
            result = -V_od[inner_idx]
        else
            # Outer generators = first K-1 indices, innermost = -V_od[last index]
            outer_indices = gen_is[1:(end-1)]
            outer_generators = QuExpr[S[i] for i in outer_indices]
            result = multi_nested_commutator(outer_generators, -V_od[inner_idx])
        end
    end

    # Cache result
    lock(_bch_cache_lock) do
        _bch_cache[cache_key] = result
    end
    return result
end

"""
    _ordered_compositions(target::Int, k::Int, max_val::Int, S::Vector{QuExpr})

Generate all ordered k-tuples of positive integers summing to target,
where each value is ≤ max_val and appears in S (is assigned).

Returns a materialized vector of index tuples.
"""
function _ordered_compositions(target::Int, k::Int, max_val::Int, S::Vector{QuExpr})
    result = Vector{Vector{Int}}()
    _ordered_compositions_rec!(result, Int[], target, k, max_val, S)
    return result
end

function _ordered_compositions_rec!(
    result::Vector{Vector{Int}},
    current::Vector{Int},
    remaining::Int,
    remaining_slots::Int,
    max_val::Int,
    S::Vector{QuExpr},
)
    if remaining_slots == 0
        if remaining == 0
            push!(result, copy(current))
        end
        return
    end

    # Minimum value needed: at least remaining - (remaining_slots - 1) * max_val
    # Maximum value allowed: min(max_val, remaining - (remaining_slots - 1) * 1)
    min_val = max(1, remaining - (remaining_slots - 1) * max_val)
    max_avail = min(max_val, remaining - (remaining_slots - 1))

    for val = min_val:max_avail
        if isassigned(S, val)
            push!(current, val)
            _ordered_compositions_rec!(
                result, current, remaining - val, remaining_slots - 1, max_val, S
            )
            pop!(current)
        end
    end
end

"""
    _collect_bch_terms_at_order_sequential(S, V_od, n, H0, V, P; parallel=false)

Collect all BCH expansion terms at order n using sequential, memoized computation.

Unlike the original compositions-based approach, this function:
1. Enumerates only the UNIQUE generator-indexed nested commutators needed
2. Caches all intermediate results (e.g., [S₁, V] computed once, reused in [S₂, [S₁, V]])
3. Avoids re-computing the same commutator recursively

The BCH expansion at order n involves:
- For each k ≥ 1 nested commutators, generator orders i₁+...+iₖ summing to n-1 (V-base) or n (H0-base)
- BCH coefficient 1/k!
"""
function _collect_bch_terms_at_order_sequential(
    S::Vector{QuExpr},
    V_od::Vector{QuExpr},
    n::Int,
    H0::QuExpr,
    V::QuExpr,
    P::Subspace;
    parallel::Bool = false,
)
    result = QuExpr()
    max_gen = n - 1  # only S₁...S_{n-1} are available

    # Collect all unique generator-index tuples needed
    # Structure: (indices_tuple, base_type) pairs
    # base_type is :V or ("H0", inner_idx)

    work_specs = Vector{Tuple{Vector{Int},Union{Symbol,Tuple{Symbol,Int}}}}()

    # Precompute BCH coefficients 1/k! for all k needed
    coeffs = Dict{Int,Rational{BigInt}}(k => big(1) // factorial(big(k)) for k = 1:n)

    # Enumerate all k (nesting depth) and all permutations of generator indices
    # that sum to the correct target
    for k = 1:n
        # V-based contributions: Σ i_j = n - 1, each i_j ≤ max_gen
        target_V = n - 1
        if target_V >= k
            # Generate all ordered k-tuples of positive integers summing to target_V
            for indices in _ordered_compositions(target_V, k, max_gen, S)
                push!(work_specs, (indices, :V))
            end
        end

        # H0-based contributions: Σ i_j = n, each i_j ≤ max_gen
        # The innermost generator index determines the V_od term used
        target_H0 = n
        if target_H0 >= k
            for indices in _ordered_compositions(target_H0, k, max_gen, S)
                inner_idx = indices[end]
                if isassigned(V_od, inner_idx) && !isempty(V_od[inner_idx].terms)
                    push!(work_specs, (indices, (:H0, inner_idx)))
                end
            end
        end
    end

    # Compute each unique contribution with memoization
    num_specs = length(work_specs)
    if parallel && n > 3 && Threads.nthreads() > 1 && num_specs > 1
        results = Vector{QuExpr}(undef, num_specs)
        Threads.@threads for i = 1:num_specs
            idx, base_type = work_specs[i]
            k = length(idx)
            base_key = base_type isa Symbol ? "V" : "H0:$(base_type[2])"
            results[i] = coeffs[k] * _nested_comm_memoized(
                Tuple(idx), base_key, S, V, V_od
            )
        end
        # Accumulate without normal_form to avoid Symbolics explosion
        result_terms = Dict{QuTerm,Number}()
        for r in results
            for (term, coeff) in r.terms
                result_terms[term] = get(result_terms, term, 0) + coeff
            end
        end
        result = QuExpr(result_terms)
    else
        result_terms = Dict{QuTerm,Number}()
        for (idx, base_type) in work_specs
            k = length(idx)
            base_key = base_type isa Symbol ? "V" : "H0:$(base_type[2])"
            term = coeffs[k] * _nested_comm_memoized(
                Tuple(idx), base_key, S, V, V_od
            )
            for (t, c) in term.terms
                result_terms[t] = get(result_terms, t, 0) + c
            end
        end
        result = QuExpr(result_terms)
    end

    if is_fastmode_flat()
        result = _collapse_transition_terms(result)
    end

    return result
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
    sw_generator(H::QuExpr, P::Subspace; order::Int=1, kwargs...)

Compute only the generator S for the Schrieffer-Wolff transformation,
without computing the full effective Hamiltonian.

This is a convenience wrapper around `schrieffer_wolff` that returns only
the generator S. It delegates to the full BCH-based algorithm, so it
produces correct generators at any order.

# Arguments
- `H`: The full Hamiltonian
- `P`: The subspace defining the block-diagonal structure
- `order`: Perturbation order for S (default: 1)
- `include_QQ`: Whether to treat Q↔Q couplings as perturbations (default: true)
- `simplify_mode`: Simplification mode for final output (default: `:fast`)
- `simplify_generator`: Whether to simplify the generator S (default: false)
- `parallel`: Use multi-threading for orders > 3 (default: false)

# Returns
- The generator S of the transformation
"""
function sw_generator(
    H::QuExpr,
    P::Subspace;
    order::Int = 1,
    include_QQ::Bool = true,
    simplify_mode::Symbol = :fast,
    simplify_generator::Bool = false,
    parallel::Bool = false,
    fastmode_flat::Bool = false,
)
    H = normal_form(H)
    H_d, H_od = decompose(H, P)

    # Extract Q↔Q couplings if requested
    H0 = H_d
    V_total = H_od
    if include_QQ
        V_QQ = _extract_QQ_terms(H_d, P)
        if !isempty(V_QQ.terms)
            H0 = normal_form(H_d - V_QQ)
            V_total = normal_form(H_od + V_QQ)
        end
    end

    # For order 1, compute S₁ directly (the iterative approach works for n=1)
    if order == 1
        _, V_od1 = decompose(V_total, P)
        S1 = solve_for_generator(H0, V_od1, P)
        if simplify_generator
            return simplify_coefficients(S1; mode = simplify_mode)
        else
            return S1
        end
    end

    # For order ≥ 2, delegate to the full BCH-based algorithm
    result = schrieffer_wolff(
        H, P;
        order = order,
        diagonal_only = false,
        include_QQ = include_QQ,
        simplify_mode = simplify_mode,
        simplify_generator = simplify_generator,
        parallel = parallel,
        fastmode_flat = fastmode_flat,
    )
    return result.S
end

"""
    project_to_subspace(H::QuExpr, P::Subspace)

Project an operator onto the subspace P.

This replaces diagonal operators by their eigenvalues in P:
- σz → eigenvalue (e.g., -1 for spin down)
- σ⁺σ⁻ → 0 for spin down, 1 for spin up
- a†a → eigenvalue (e.g., 0 for vacuum)

And removes any remaining off-diagonal terms.

Returns a projected `QuExpr`. For performance, the return value is not
automatically passed through `normal_form`; callers that require canonical
ordering should apply `normal_form` explicitly.
"""
function project_to_subspace(H::QuExpr, P::Subspace)
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

        # Check if this is a transition operator constraint (N-level system projector)
        trans_info = get_transition_constraint_info(constraint)
        if trans_info !== nothing
            result = substitute_transition_projections(result, trans_info)
            continue
        end

        # Check if this is a number operator constraint (a†a => n)
        if is_number_constraint(constraint)
            result = substitute_number_operator_projection(result, constraint)
            continue
        end

        # For other constraints, try direct substitution
        result = substitute_operator_eigenvalue(
            result,
            constraint.operator,
            constraint.eigenvalue,
        )
    end

    return result
end

"""
    substitute_number_operator_projection(H::QuExpr, constraint::OperatorConstraint)

Substitute boson/fermion operators for a number operator constraint a†a => n.

For n = 0 (vacuum):
- a†a → 0
- Any term containing a or a† vanishes (they take us out of vacuum)

For n > 0:
- a†a → n
- Terms with unbalanced a†/a are removed by diagonal_part already
"""
function substitute_number_operator_projection(H::QuExpr, constraint::OperatorConstraint)
    n = constraint.eigenvalue

    # Extract the operator info from the constraint
    op_term = first(constraint.operator.terms)[1]
    ops = op_term.bares.v
    # For a†a, ops[1] is creation, ops[2] is annihilation
    op_name = ops[1].name
    op_inds = ops[1].inds
    is_boson = ops[1].t == BosonCreate_
    create_type = is_boson ? BosonCreate_ : FermionCreate_
    destroy_type = is_boson ? BosonDestroy_ : FermionDestroy_

    result_terms = Dict{QuTerm,Number}()

    for (term, coeff) in H.terms
        term_ops = term.bares.v
        new_ops = BaseOperator[]
        eigenvalue_factor = 1
        term_vanishes = false
        i = 1

        while i <= length(term_ops)
            op = term_ops[i]

            # Check if this operator matches our constraint
            if op.name == op_name && op.inds == op_inds
                if op.t == create_type || op.t == destroy_type
                    if n == 0
                        # In vacuum, any a or a† kills the term
                        term_vanishes = true
                        break
                    else
                        # For n > 0, keep the operator (balanced terms survive)
                        push!(new_ops, op)
                    end
                elseif i < length(term_ops)
                    # Check for a†a pattern
                    next_op = term_ops[i+1]
                    if op.t == create_type &&
                       next_op.t == destroy_type &&
                       next_op.name == op_name &&
                       next_op.inds == op_inds
                        # Found a†a - replace with eigenvalue n
                        eigenvalue_factor *= n
                        i += 2
                        continue
                    else
                        push!(new_ops, op)
                    end
                else
                    push!(new_ops, op)
                end
            else
                push!(new_ops, op)
            end
            i += 1
        end

        if term_vanishes
            continue
        end

        # Build the new term
        if isempty(new_ops)
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

        # If eigenvalue factor is 0, skip the term
        if eigenvalue_factor == 0
            continue
        end

        new_coeff = coeff * eigenvalue_factor
        result_terms[new_term] = get(result_terms, new_term, 0) + new_coeff
    end

    if is_fastmode_flat()
        return QuExpr(result_terms)
    end
    return normal_form(QuExpr(result_terms))
end

"""
    substitute_spin_projection(H::QuExpr, spin_name, spin_inds, pm_eigenvalue)

Substitute σ⁺σ⁻ (or σ⁺(i)σ⁻(i)) by its eigenvalue in the spin subspace.
pm_eigenvalue is 0 for spin-down, 1 for spin-up.
"""
function substitute_spin_projection(H::QuExpr, spin_name, spin_inds, pm_eigenvalue)
    result_terms = Dict{QuTerm,Number}()

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

        new_coeff = coeff * eigenvalue_factor
        result_terms[new_term] = get(result_terms, new_term, 0) + new_coeff
    end

    if is_fastmode_flat()
        return QuExpr(result_terms)
    end
    return normal_form(QuExpr(result_terms))
end

"""
    substitute_transition_projections(H::QuExpr, trans_info::NamedTuple)

Substitute transition operators by their expectation values for N-level systems.

For a constraint on projector |k⟩⟨k| with eigenvalue 1 (we ARE in state k):
- |k⟩⟨k| → 1 (the constrained projector)
- |m⟩⟨m| for m ≠ k → 0 (orthogonal projectors)
- |i⟩⟨j| where neither i nor j is k → 0 (off-diagonal, ⟨k|i⟩⟨j|k⟩ = 0)
- |k⟩⟨m| or |m⟩⟨k| for m ≠ k → 0 (these are off-diagonal in the projected subspace)

For eigenvalue 0 (we are NOT in state k):
- |k⟩⟨k| → 0

trans_info contains: name, inds, N, state (the constrained state k), eigenvalue
"""
function substitute_transition_projections(H::QuExpr, trans_info::NamedTuple)
    result_terms = Dict{QuTerm,Number}()

    constrained_state = trans_info.state
    eigenvalue = trans_info.eigenvalue
    op_name = trans_info.name
    op_inds = trans_info.inds
    N = trans_info.N

    for (term, coeff) in H.terms
        ops = term.bares.v
        new_ops = BaseOperator[]
        eigenvalue_factor = 1
        term_vanishes = false

        for op in ops
            # Check if this is a transition operator with matching name/inds
            if op.t == Transition_ && op.name == op_name && op.inds == op_inds
                indices = get_transition_indices(op)
                if indices !== nothing
                    i, j = indices

                    if eigenvalue == 1
                        # We ARE in state `constrained_state` (call it k)
                        # Projection: P = |k⟩⟨k|
                        # P |i⟩⟨j| P = |k⟩⟨k|i⟩⟨j|k⟩⟨k| = δ_{ki} δ_{jk} |k⟩⟨k|
                        # So only |k⟩⟨k| survives, and it becomes the identity in the subspace

                        if i == constrained_state && j == constrained_state
                            # |k⟩⟨k| → 1 (identity in the projected subspace)
                            eigenvalue_factor *= 1
                            # Don't push the operator - it's replaced by 1
                        else
                            # Any other transition operator gives 0
                            term_vanishes = true
                            break
                        end
                    else  # eigenvalue == 0
                        # We are NOT in state `constrained_state`
                        if i == constrained_state && j == constrained_state
                            # |k⟩⟨k| → 0
                            term_vanishes = true
                            break
                        else
                            # Other operators stay (we don't know exactly which state we're in)
                            push!(new_ops, op)
                        end
                    end
                    continue
                end
            end
            # Non-matching operator: keep it
            push!(new_ops, op)
        end

        # Skip term if it vanishes
        if term_vanishes
            continue
        end

        # Build the new term
        if isempty(new_ops)
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

        new_coeff = coeff * eigenvalue_factor
        result_terms[new_term] = get(result_terms, new_term, 0) + new_coeff
    end

    if is_fastmode_flat()
        return QuExpr(result_terms)
    end
    return normal_form(QuExpr(result_terms))
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

    result_terms = Dict{QuTerm,Number}()

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
                new_coeff = coeff * eigenvalue_factor
                result_terms[new_term] = get(result_terms, new_term, 0) + new_coeff
                continue
            end
        end

        if !matched
            result_terms[term] = get(result_terms, term, 0) + coeff
        end
    end

    if is_fastmode_flat()
        return QuExpr(result_terms)
    end
    return normal_form(QuExpr(result_terms))
end
