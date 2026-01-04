"""
    Baker-Campbell-Hausdorff (BCH) expansion for unitary transformations.

Computes e^S H e^{-S} = H + [S,H] + (1/2!)[S,[S,H]] + (1/3!)[S,[S,[S,H]]] + ...
"""

export commutator_series, nested_commutator, bch_transform

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
    for _ in 1:n
        result = comm(S, result)
        result = normal_form(result)
    end
    return result
end

"""
    commutator_series(S::QuExpr, H::QuExpr, order::Int)

Compute the BCH expansion of e^S H e^{-S} to the given order.

e^S H e^{-S} = Σₙ (1/n!) [S, [S, [..., [S, H]...]]]  (n nested commutators)

The expansion is truncated at `order` nested commutators.

# Arguments
- `S`: The generator of the unitary transformation (anti-Hermitian: S† = -S)
- `H`: The Hamiltonian to transform
- `order`: Maximum number of nested commutators to include

# Returns
- The transformed Hamiltonian as a QuExpr
"""
function commutator_series(S::QuExpr, H::QuExpr, order::Int)
    order >= 0 || throw(ArgumentError("order must be non-negative, got $order"))
    
    result = QuExpr()
    factorial_n = 1
    
    for n in 0:order
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
    bch_transform(S::QuExpr, H::QuExpr; order::Int=4)

Compute e^S H e^{-S} using the BCH formula to the specified order.

This is an alias for `commutator_series` with a more intuitive name.
"""
bch_transform(S::QuExpr, H::QuExpr; order::Int=4) = commutator_series(S, H, order)
