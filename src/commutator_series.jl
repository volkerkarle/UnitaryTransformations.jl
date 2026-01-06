"""
    Baker-Campbell-Hausdorff (BCH) expansion for unitary transformations.

Provides several BCH-related functions:
- `bch_transform(S, A; order)`: Compute e^S A e^{-S} (adjoint action)
- `bch_combine(A, B; order)`: Compute Z where e^A e^B = e^Z
- `nested_commutator(S, H, n)`: Compute n-fold nested commutator
"""

export commutator_series,
    nested_commutator, multi_nested_commutator, compositions, bch_transform, bch_combine

using QuantumAlgebra
using QuantumAlgebra: QuExpr, normal_form, comm

"""
    nested_commutator(S::QuExpr, H::QuExpr, n::Int)

Compute the n-fold nested commutator [S, [S, [..., [S, H]...]]]
with S appearing n times.

- n=0 returns H
- n=1 returns [S, H]
- n=2 returns [S, [S, H]]
- etc.
"""
function nested_commutator(S::QuExpr, H::QuExpr, n::Int)
    n >= 0 || throw(ArgumentError("n must be non-negative, got $n"))

    result = H
    for _ = 1:n
        result = comm(S, result)
        result = normal_form(result)
    end
    return result
end

"""
    multi_nested_commutator(generators::Vector{QuExpr}, X::QuExpr)

Compute the nested commutator [g₁, [g₂, [..., [gₖ, X]...]]] for a list of generators.

The generators are applied from right to left (innermost first), so:
- `multi_nested_commutator([S₁, S₂], X)` computes `[S₁, [S₂, X]]`
- `multi_nested_commutator([S₁], X)` computes `[S₁, X]`
- `multi_nested_commutator(QuExpr[], X)` returns `X`

# Arguments
- `generators`: Vector of generator expressions (can have different orders)
- `X`: The base operator (innermost argument)

# Returns
- The nested commutator as a QuExpr
"""
function multi_nested_commutator(generators::Vector{QuExpr}, X::QuExpr)
    isempty(generators) && return X

    result = X
    # Apply from right to left: for [S₁, [S₂, X]], we first compute [S₂, X], then [S₁, ...]
    for g in reverse(generators)
        result = normal_form(comm(g, result))
    end
    return result
end

"""
    compositions(n::Int, k::Int; min_val::Int=1, max_val::Int=n)

Generate all k-tuples of positive integers that sum to n.

A composition is an ordered partition: (1,2) and (2,1) are different compositions of 3 into 2 parts.

# Arguments
- `n`: The target sum
- `k`: The number of parts (length of each tuple)
- `min_val`: Minimum value for each part (default: 1)
- `max_val`: Maximum value for each part (default: n)

# Returns
- Vector of Vector{Int}, each inner vector is a composition

# Examples
```julia
compositions(3, 2)  # [(1,2), (2,1)]
compositions(4, 2)  # [(1,3), (2,2), (3,1)]
compositions(4, 2; max_val=2)  # [(2,2)]
compositions(4, 3)  # [(1,1,2), (1,2,1), (2,1,1)]
```
"""
function compositions(n::Int, k::Int; min_val::Int = 1, max_val::Int = n)
    n >= 0 || throw(ArgumentError("n must be non-negative, got $n"))
    k >= 0 || throw(ArgumentError("k must be non-negative, got $k"))

    if k == 0
        return n == 0 ? [Int[]] : Vector{Int}[]
    end

    result = Vector{Vector{Int}}()

    function generate(current::Vector{Int}, remaining::Int, parts_left::Int)
        if parts_left == 0
            if remaining == 0
                push!(result, copy(current))
            end
            return
        end

        # For the current part, try all valid values
        # Must leave room for remaining parts to have at least min_val each
        min_remaining = min_val * (parts_left - 1)
        max_remaining = max_val * (parts_left - 1)

        lo = max(min_val, remaining - max_remaining)
        hi = min(max_val, remaining - min_remaining)

        for val = lo:hi
            push!(current, val)
            generate(current, remaining - val, parts_left - 1)
            pop!(current)
        end
    end

    generate(Int[], n, k)
    return result
end

"""
    commutator_series(S::QuExpr, H::QuExpr, order::Int)

Compute the BCH expansion of e^S H e^{-S} to the given order.

e^S H e^{-S} = Σₙ (1/n!) [S, [S, [..., [S, H]...]]]  (n nested commutators)

The expansion is truncated at `order` nested commutators.

# Arguments
- `S`: The generator of the unitary transformation (anti-Hermitian: S† = -S)
- `H`: The operator to transform
- `order`: Maximum number of nested commutators to include

# Returns
- The transformed operator as a QuExpr
"""
function commutator_series(S::QuExpr, H::QuExpr, order::Int)
    order >= 0 || throw(ArgumentError("order must be non-negative, got $order"))

    result = QuExpr()
    factorial_n = 1

    for n = 0:order
        if n > 0
            factorial_n *= n
        end

        # Compute [S, [S, [..., [S, H]...]]] with n nested commutators
        term = nested_commutator(S, H, n)

        # Add (1/n!) * term to result
        result = result + term * (1 // factorial_n)
        result = normal_form(result)
    end

    return result
end

"""
    bch_transform(S::QuExpr, A::QuExpr; order::Int=4)

Compute e^S A e^{-S} using the BCH formula to the specified order.

This computes the adjoint action of e^S on an operator A:
    e^S A e^{-S} = A + [S,A] + (1/2!)[S,[S,A]] + (1/3!)[S,[S,[S,A]]] + ...

# Arguments
- `S`: The generator (typically anti-Hermitian: S† = -S)
- `A`: The operator to transform (e.g., a Hamiltonian)
- `order`: Number of nested commutators to include (default: 4)

# Example
```julia
using QuantumAlgebra, UnitaryTransformations, Symbolics
@variables g Δ

# Generator from SW transformation
S = (g/Δ) * (a'()*σm() - a()*σp())

# Transform the Hamiltonian
H = Δ/2 * σz() + g * (a'()*σm() + a()*σp())
H_transformed = bch_transform(S, H; order=2)
```
"""
bch_transform(S::QuExpr, A::QuExpr; order::Int = 4) = commutator_series(S, A, order)

"""
    bch_combine(A::QuExpr, B::QuExpr; order::Int=4)

Compute Z such that e^A e^B = e^Z using the BCH formula.

The BCH formula gives:
    Z = A + B + (1/2)[A,B] + (1/12)[A,[A,B]] - (1/12)[B,[A,B]] 
        + (1/24)[B,[A,[A,B]]] + ...

This is useful for combining sequential unitary transformations.

# Arguments
- `A`: First generator
- `B`: Second generator  
- `order`: Order of the expansion (default: 4)
  - order=1: Z = A + B
  - order=2: Z = A + B + (1/2)[A,B]
  - order=3: Z = A + B + (1/2)[A,B] + (1/12)[A,[A,B]] - (1/12)[B,[A,B]]
  - etc.

# Returns
- Z such that e^A e^B = e^Z (to the specified order)

# Example
```julia
using QuantumAlgebra, UnitaryTransformations

# Combine two rotations
A = param(:θ) * σx()
B = param(:φ) * σy()
Z = bch_combine(A, B; order=3)
# e^A e^B ≈ e^Z
```

# Notes
The full BCH series involves Bernoulli numbers and can be written as:
    Z = A + B + Σₙ (Bₙ/n!) Cₙ(A,B)
where Cₙ are nested commutators and Bₙ are Bernoulli numbers.

For practical purposes, we compute terms explicitly up to the given order.
"""
function bch_combine(A::QuExpr, B::QuExpr; order::Int = 4)
    order >= 1 || throw(ArgumentError("order must be at least 1, got $order"))

    # Order 1: Z = A + B
    Z = normal_form(A + B)

    if order == 1
        return Z
    end

    # Order 2: Z += (1/2)[A,B]
    C_AB = normal_form(comm(A, B))  # [A,B]
    Z = normal_form(Z + C_AB * (1 // 2))

    if order == 2
        return Z
    end

    # Order 3: Z += (1/12)[A,[A,B]] - (1/12)[B,[A,B]]
    #        = (1/12)([A,[A,B]] + [[A,B],B])
    C_A_AB = normal_form(comm(A, C_AB))  # [A,[A,B]]
    C_AB_B = normal_form(comm(C_AB, B))  # [[A,B],B]
    Z = normal_form(Z + (C_A_AB + C_AB_B) * (1 // 12))

    if order == 3
        return Z
    end

    # Order 4: Z += -(1/24)[B,[A,[A,B]]]
    #        Note: [A,[B,[A,B]]] term vanishes by Jacobi identity rearrangement
    C_B_A_AB = normal_form(comm(B, C_A_AB))  # [B,[A,[A,B]]]
    Z = normal_form(Z - C_B_A_AB * (1 // 24))

    if order == 4
        return Z
    end

    # Order 5 and beyond: use recursive formula
    # The general BCH formula involves Bernoulli numbers and gets complex.
    # For higher orders, we use the Dynkin formula approach.

    if order >= 5
        # Higher order terms using the recursive structure
        # [A,[A,[A,[A,B]]]] terms, [B,[B,[A,B]]] terms, mixed terms...

        # Order 5 terms:
        # -(1/720)[A,[A,[A,[A,B]]]] - (1/720)[B,[B,[B,[A,B]]]]
        # +(1/360)[A,[B,[B,[A,B]]]] + (1/360)[B,[A,[A,[A,B]]]]
        # +(1/120)[A,[B,[A,[A,B]]]] + (1/120)[B,[A,[B,[A,B]]]]

        C_A_A_AB = normal_form(comm(A, C_A_AB))  # [A,[A,[A,B]]]
        C_B_AB = normal_form(comm(B, C_AB))  # [B,[A,B]]
        C_B_B_AB = normal_form(comm(B, C_B_AB))  # [B,[B,[A,B]]]

        # [A,[A,[A,[A,B]]]]
        C_A_A_A_AB = normal_form(comm(A, C_A_A_AB))
        # [B,[B,[B,[A,B]]]]
        C_B_B_B_AB = normal_form(comm(B, C_B_B_AB))

        Z = normal_form(Z - (C_A_A_A_AB + C_B_B_B_AB) * (1 // 720))

        # Additional order 5 cross terms
        C_A_B_B_AB = normal_form(comm(A, C_B_B_AB))
        C_B_A_A_AB = normal_form(comm(B, C_A_A_AB))
        Z = normal_form(Z + (C_A_B_B_AB + C_B_A_A_AB) * (1 // 360))

        C_A_B_A_AB = normal_form(comm(A, C_B_A_AB))
        C_B_A_B_AB = normal_form(comm(B, normal_form(comm(A, C_B_AB))))
        Z = normal_form(Z + (C_A_B_A_AB + C_B_A_B_AB) * (1 // 120))
    end

    if order >= 6
        # For order 6+, the formulas get very complex.
        # A general implementation would use the Goldberg coefficients
        # or compute via the integral formula.
        @warn "BCH order > 5 uses approximate coefficients; results may be incomplete"
    end

    return Z
end
