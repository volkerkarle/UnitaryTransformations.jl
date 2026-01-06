"""
    Magnus Expansion for periodically driven quantum systems.

The Magnus expansion provides a way to compute an effective time-independent
Hamiltonian for a periodically driven system H(t) = H(t + T).

For a time-dependent Schrödinger equation:
    i ∂U/∂t = H(t) U(t)

The solution can be written as U(t) = exp(Ω(t)) where Ω is the Magnus series:
    Ω(t) = Ω₁(t) + Ω₂(t) + Ω₃(t) + ...

For periodic driving, the effective (Floquet) Hamiltonian after one period is:
    H_eff = i Ω(T) / T

This module handles Hamiltonians in Fourier representation:
    H(t) = Σₙ Hₙ e^{inωt}

where ω = 2π/T is the driving frequency.

The implementation supports arbitrary order expansions using a recursive
algorithm based on the high-frequency expansion formula.
"""

export magnus_expansion, check_hermiticity, FourierHamiltonian
export MagnusResult

using QuantumAlgebra
using QuantumAlgebra: QuExpr, normal_form, comm

using Symbolics
using Symbolics: Num

import ..UnitaryTransformations: simplify_coefficients

# =============================================================================
# FourierHamiltonian type
# =============================================================================

"""
    FourierHamiltonian

Represents a time-periodic Hamiltonian in Fourier representation:
    H(t) = Σₙ Hₙ e^{inωt}

The modes are stored as a Dict{Int, QuExpr} mapping Fourier index n to operator Hₙ.

For Hermiticity, we require H₋ₙ = Hₙ†.
"""
struct FourierHamiltonian
    modes::Dict{Int,QuExpr}
    ω::Num  # Driving frequency (symbolic)

    function FourierHamiltonian(modes::Dict{Int,QuExpr}, ω)
        ω_num = ω isa Num ? ω : Num(ω)
        new(modes, ω_num)
    end
end

"""
    FourierHamiltonian(modes::Dict, ω)

Construct a Fourier Hamiltonian from a dictionary of modes and frequency.

# Example
```julia
@variables Δ Ω ω
modes = Dict(
    0  => Δ/2 * σz(),
    1  => Ω/4 * (σp() + σm()),
    -1 => Ω/4 * (σp() + σm())
)
H = FourierHamiltonian(modes, ω)
```
"""
FourierHamiltonian(modes::Dict, ω) = FourierHamiltonian(convert(Dict{Int,QuExpr}, modes), ω)

Base.getindex(H::FourierHamiltonian, n::Int) = get(H.modes, n, zero(QuExpr))
Base.haskey(H::FourierHamiltonian, n::Int) = haskey(H.modes, n)
Base.keys(H::FourierHamiltonian) = keys(H.modes)

# =============================================================================
# Hermiticity check
# =============================================================================

"""
    check_hermiticity(H::FourierHamiltonian; warn_only::Bool=false)

Check that the Fourier Hamiltonian satisfies H₋ₙ = Hₙ† for all modes.

Returns `true` if Hermitian, throws an error otherwise.
"""
function check_hermiticity(H::FourierHamiltonian; warn_only::Bool = false)
    non_hermitian_modes = Int[]

    for n in keys(H.modes)
        if n == 0
            H0 = H.modes[0]
            H0_dag = normal_form(H0')
            if normal_form(H0) != H0_dag
                push!(non_hermitian_modes, 0)
            end
        elseif n > 0
            Hn = H.modes[n]
            Hn_dag = normal_form(Hn')
            H_minus_n = get(H.modes, -n, zero(QuExpr))

            if normal_form(H_minus_n) != Hn_dag
                push!(non_hermitian_modes, n)
            end
        end
    end

    if !isempty(non_hermitian_modes)
        msg =
            "Fourier Hamiltonian is not Hermitian. " *
            "Modes violating H₋ₙ = Hₙ†: $non_hermitian_modes"
        if warn_only
            @warn msg
            return false
        else
            throw(ArgumentError(msg))
        end
    end

    return true
end

function check_hermiticity(modes::Dict{Int,QuExpr}; warn_only::Bool = false)
    @variables _ω_placeholder
    H = FourierHamiltonian(modes, _ω_placeholder)
    return check_hermiticity(H; warn_only = warn_only)
end

# =============================================================================
# Nested commutator computation
# =============================================================================

"""
    nested_commutator(operators::Vector{QuExpr})

Compute the left-nested commutator [...[[A₁, A₂], A₃], ..., Aₙ].
"""
function nested_commutator(operators::Vector{QuExpr})
    isempty(operators) && return zero(QuExpr)
    length(operators) == 1 && return operators[1]

    result = operators[1]
    for i = 2:length(operators)
        result = normal_form(comm(result, operators[i]))
    end
    return result
end

"""
    nested_commutator(H::FourierHamiltonian, indices)

Compute the nested commutator for given Fourier indices.
"""
function nested_commutator(H::FourierHamiltonian, indices)
    for n in indices
        if !haskey(H.modes, n)
            return nothing
        end
    end
    operators = [H.modes[n] for n in indices]
    return nested_commutator(operators)
end

# =============================================================================
# Dynamic Magnus expansion - Recursive formula
# =============================================================================

"""
    _magnus_order_k(H::FourierHamiltonian, k::Int; cache=nothing)

Compute the k-th order Magnus term dynamically using the recursive formula.

The high-frequency/Floquet-Magnus expansion gives:

    H_eff^(1) = H₀

    H_eff^(k) = Σ (over valid index sets) C(n₁,...,nₖ) × [...[[Hₙ₁, Hₙ₂], Hₙ₃],..., Hₙₖ]

where the sum is over indices satisfying:
- n₁ + n₂ + ... + nₖ = 0 (resonance condition)
- The combination is not reducible to a lower order

The coefficients C(n₁,...,nₖ) are computed from nested time integrals.

For practical computation, we use the explicit formulas:
- Order 2: Σₙ>₀ -[Hₙ, H₋ₙ]/(nω)
- Order k≥3: Recursive structure based on partial sums
"""
function _magnus_order_k(H::FourierHamiltonian, k::Int)
    k >= 1 || throw(ArgumentError("Order must be at least 1"))

    if k == 1
        return H[0]
    end

    if k == 2
        return _magnus_order_2(H)
    end

    # For k ≥ 3, use the general recursive approach
    return _magnus_higher_order(H, k)
end

"""
    _magnus_order_2(H::FourierHamiltonian)

Second-order Magnus: H_eff^(2) = Σₙ>₀ -[Hₙ, H₋ₙ]/(nω)
"""
function _magnus_order_2(H::FourierHamiltonian)
    result = zero(QuExpr)
    ω = H.ω

    for n in keys(H.modes)
        if n > 0 && haskey(H.modes, -n)
            Hn = H.modes[n]
            Hmn = H.modes[-n]
            comm_term = normal_form(comm(Hn, Hmn))
            coeff = -1 / (n * ω)
            result = normal_form(result + coeff * comm_term)
        end
    end

    return result
end

"""
    _magnus_higher_order(H::FourierHamiltonian, k::Int)

Compute Magnus term for order k ≥ 3 using systematic enumeration.

The effective Hamiltonian at order k involves nested commutators of k Fourier
components with indices (n₁, ..., nₖ) satisfying:
1. n₁ + n₂ + ... + nₖ = 0 (resonance)
2. Not reducible to lower order (no adjacent identical non-zero indices)

The coefficient for each term depends on the partial sums of indices.
"""
function _magnus_higher_order(H::FourierHamiltonian, k::Int)
    result = zero(QuExpr)
    ω = H.ω

    mode_indices = sort(collect(keys(H.modes)))
    positive_indices = filter(n -> n > 0, mode_indices)

    # Generate ALL valid index combinations (including those with zeros)
    all_combinations = _generate_all_valid_combinations(mode_indices, k)

    for indices in all_combinations
        # Check canonical ordering (first non-zero should be positive)
        if !_is_canonical_ordering(indices, positive_indices)
            continue
        end

        # Skip if adjacent identical non-zero indices (commutator vanishes)
        if _has_adjacent_identical_nonzero(indices)
            continue
        end

        # Compute coefficient
        coeff = _compute_coefficient_general(indices, ω)
        if coeff === nothing
            continue
        end

        # Compute nested commutator
        comm_result = nested_commutator(H, indices)
        if comm_result === nothing || comm_result == zero(QuExpr)
            continue
        end

        result = normal_form(result + coeff * comm_result)
    end

    return result
end

"""
    _generate_all_valid_combinations(mode_indices, k)

Generate all k-tuples of indices that sum to zero.
"""
function _generate_all_valid_combinations(mode_indices, k)
    result = Vector{Vector{Int}}()
    _gen_all_combinations!(result, mode_indices, k, Int[], 0)
    return result
end

function _gen_all_combinations!(result, indices, remaining, current, current_sum)
    if remaining == 1
        needed = -current_sum
        if needed in indices
            push!(result, vcat(current, [needed]))
        end
        return
    end
    for idx in indices
        push!(current, idx)
        _gen_all_combinations!(result, indices, remaining - 1, current, current_sum + idx)
        pop!(current)
    end
end

"""
Check if indices contain adjacent identical non-zero values.
"""
function _has_adjacent_identical_nonzero(indices)
    for i = 1:(length(indices)-1)
        if indices[i] == indices[i+1] && indices[i] != 0
            return true
        end
    end
    return false
end

"""
    _compute_coefficient_general(indices, ω)

Compute coefficient for arbitrary index combination using the general formula.

For indices (n₁, ..., nₖ) with Σnᵢ = 0, the coefficient is:
    C = 1 / (ω^(k-1) × Π_{j=1}^{k-1} sⱼ)

where sⱼ = n₁ + ... + nⱼ are partial sums.

Returns `nothing` if the term is reducible (partial sum is zero at intermediate position).
"""
function _compute_coefficient_general(indices, ω)
    k = length(indices)

    # Compute partial sums
    partial_sums = cumsum(indices)

    # Check for zeros in partial sums (excluding the last which is always zero)
    # A zero at position j means n₁ + ... + nⱼ = 0, which makes this a reducible term
    # (it factorizes into lower-order contributions)
    for j = 1:(k-1)
        if partial_sums[j] == 0
            return nothing  # Reducible term
        end
    end

    # All partial sums are non-zero - compute coefficient
    product = prod(partial_sums[1:(k-1)])

    return 1 / (ω^(k-1) * product)
end

"""
Check if an index sequence is in canonical order.

For Hermitian systems, we define canonical ordering as:
- The first non-zero index is positive
- This avoids counting both (n,...,-n) and (-n,...,n) which give equivalent results
"""
function _is_canonical_ordering(indices, positive_indices)
    for n in indices
        if n != 0
            return n > 0
        end
    end
    return true  # All zeros (shouldn't happen for our use case)
end

# =============================================================================
# Result type
# =============================================================================

"""
    MagnusResult

Result type for Magnus expansion containing H_eff and individual order contributions.
"""
struct MagnusResult
    H_eff::QuExpr
    orders::Dict{Int,QuExpr}
end

Base.getindex(r::MagnusResult, k::Int) = get(r.orders, k, zero(QuExpr))
Base.keys(r::MagnusResult) = keys(r.orders)

function Base.getproperty(r::MagnusResult, s::Symbol)
    if s == :H_eff
        return getfield(r, :H_eff)
    elseif s == :orders
        return getfield(r, :orders)
    else
        m = match(r"^Ω(\d+)$", String(s))
        if m !== nothing
            k = parse(Int, m.captures[1])
            return get(getfield(r, :orders), k, zero(QuExpr))
        end
        error("MagnusResult has no property $s")
    end
end

Base.propertynames(r::MagnusResult) = (:H_eff, :orders, Symbol.(["Ω$k" for k = 1:10])...)

# Iteration for backward compatibility
function Base.iterate(r::MagnusResult, state = 1)
    if state == 1
        return (r.H_eff, 2)
    elseif state <= 5
        return (get(r.orders, state - 1, zero(QuExpr)), state + 1)
    else
        return nothing
    end
end
Base.length(::MagnusResult) = 5

# =============================================================================
# Main function
# =============================================================================

"""
    magnus_expansion(H::FourierHamiltonian; order::Int=2, check_hermitian::Bool=true)

Compute the Magnus expansion for a periodically driven system to arbitrary order.

# Arguments
- `H`: A `FourierHamiltonian` representing H(t) = Σₙ Hₙ e^{inωt}
- `order`: Maximum order of the expansion (default: 2, no upper limit)
- `check_hermitian`: Whether to verify H is Hermitian (default: true)

# Returns
A `MagnusResult` containing:
- `H_eff`: The effective time-independent Hamiltonian
- `orders`: Dict mapping order k to Ωₖ
- `Ω1`, `Ω2`, ...: Direct access to individual contributions

# Example
```julia
using QuantumAlgebra, UnitaryTransformations, Symbolics
QuantumAlgebra.use_σpm(true)

@variables Δ Ω ω
modes = Dict(0 => Δ/2*σz(), 1 => Ω/2*σp(), -1 => Ω/2*σm())
result = magnus_expansion(modes, ω; order=5)
println("Order 5: ", result.Ω5)
```
"""
function magnus_expansion(
    H::FourierHamiltonian;
    order::Int = 2,
    check_hermitian::Bool = true,
)
    order >= 1 || throw(ArgumentError("order must be at least 1, got $order"))

    if check_hermitian
        check_hermiticity(H)
    end

    orders = Dict{Int,QuExpr}()
    H_eff = zero(QuExpr)

    for k = 1:order
        Ωk = _magnus_order_k(H, k)
        Ωk = simplify_coefficients(Ωk)
        orders[k] = Ωk
        H_eff = normal_form(H_eff + Ωk)
    end

    H_eff = simplify_coefficients(H_eff)
    return MagnusResult(H_eff, orders)
end

function magnus_expansion(
    modes::Dict{Int,QuExpr},
    ω;
    order::Int = 2,
    check_hermitian::Bool = true,
)
    H = FourierHamiltonian(modes, ω)
    return magnus_expansion(H; order = order, check_hermitian = check_hermitian)
end

function magnus_expansion(modes::Dict{Int,T}, ω; kwargs...) where {T}
    converted = Dict{Int,QuExpr}(k => v for (k, v) in modes)
    return magnus_expansion(converted, ω; kwargs...)
end
