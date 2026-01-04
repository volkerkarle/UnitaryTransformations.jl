# [API Reference](@id api)

## Main Functions

### Schrieffer-Wolff Transformation

```@docs
schrieffer_wolff
sw_generator
```

### Subspace Definition

```@docs
Subspace
OperatorConstraint
```

### Hamiltonian Decomposition

```@docs
decompose
diagonal_part
off_diagonal_part
is_diagonal
is_off_diagonal
```

### Generator Solution

```@docs
solve_for_generator
compute_energy_denominator
```

### Projection

```@docs
project_to_subspace
```

### BCH Expansion

```@docs
bch_transform
commutator_series
nested_commutator
```

## Symbolic Utilities

### Coefficient Manipulation

```@docs
simplify_coefficients
substitute_values
extract_coefficient
collect_terms
```

### Parameter Conversion

```@docs
param_to_symbolic
clear_param_cache!
```

## Re-exported from QuantumAlgebra

The following functions are re-exported for convenience:

- `comm(A, B)` - Compute commutator ``[A, B] = AB - BA``
- `normal_form(expr)` - Normal-order an operator expression
- `a()`, `a'()` - Bosonic annihilation/creation operators
- `σx()`, `σy()`, `σz()` - Pauli matrices
- `σp()`, `σm()` - Raising/lowering operators (when `use_σpm(true)`)
- `Pr"name"` - Define symbolic parameters

## Types

### SWResult

The `schrieffer_wolff` function returns a named tuple:

```julia
result = schrieffer_wolff(H, P; order=2)

result.S      # Generator of the transformation
result.H_eff  # Block-diagonal effective Hamiltonian  
result.H_P    # H_eff projected onto subspace P
```

### Classification Enums

The package uses classification enums for operator analysis:

```julia
DIAGONAL   # Operator preserves the subspace
RAISING    # Operator raises out of subspace P
LOWERING   # Operator lowers into subspace P
MIXED      # Operator has mixed character
```
