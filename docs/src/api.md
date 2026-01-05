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

---

## Generator Solution

The core operation in Schrieffer-Wolff is solving ``[S, H_d] = -V_{od}`` for the generator ``S``.

### Main Function

```@docs
solve_for_generator
```

### Method-Specific Functions

Two methods are available, automatically selected based on the operator types:

```@docs
solve_for_generator_eigenoperator
solve_for_generator_lie
```

The **eigenoperator method** works when ``[H_d, O] = \varepsilon \cdot O`` (TLS, bosons, N-level transitions).

The **Lie algebra method** works for SU(N) systems by converting to the Cartan-Weyl basis where generators become eigenoperators.

### Supporting Functions

```@docs
compute_energy_denominator
compute_energy_eigenvalues
detect_lie_algebra_system
```

---

## Projection

```@docs
project_to_subspace
```

---

## BCH Expansion

The Baker-Campbell-Hausdorff expansion is used to compute ``e^S H e^{-S}``.

```@docs
bch_transform
commutator_series
nested_commutator
```

---

## Symbolic Utilities

### Coefficient Manipulation

```@docs
simplify_coefficients
substitute_values
extract_coefficient
collect_terms
```

### LaTeX Output

```@docs
to_latex
print_latex
show_result
```

### Parameter Conversion

Functions for converting between QuantumAlgebra's `Param` and Symbolics.jl variables:

```@docs
param_to_symbolic
symbolic_coefficient
clear_param_cache!
```

---

## Re-exported from QuantumAlgebra

The following functions are re-exported for convenience:

- `comm(A, B)` - Compute commutator ``[A, B] = AB - BA``
- `normal_form(expr)` - Normal-order an operator expression
- `a()`, `a'()` - Bosonic annihilation/creation operators
- `σx()`, `σy()`, `σz()` - Pauli matrices
- `σp()`, `σm()` - Raising/lowering operators (when `use_σpm(true)`)
- `nlevel_ops(N, name)` - N-level transition operators ``|i\rangle\langle j|``
- `su_generators(N, name)` - SU(N) generators (generalized Gell-Mann matrices)

---

## Symbolic Parameters

Use Symbolics.jl `@variables` to define symbolic parameters:

```julia
using Symbolics
@variables Δ g ω  # Define symbolic parameters
```

For N-level systems with indexed parameters:

```julia
# Create ω₁, ω₂, ..., ωₙ
ω = [Symbolics.variable(Symbol("ω", i)) for i in 1:N]
```

---

## Types

### SWResult

The `schrieffer_wolff` function returns a named tuple:

```julia
result = schrieffer_wolff(H, P; order=2)

result.S      # Generator of the transformation
result.H_eff  # Block-diagonal effective Hamiltonian  
result.H_P    # H_eff projected onto subspace P
```

### Subspace

```julia
# Single constraint
P = Subspace(σz() => -1)

# Multiple constraints
P = Subspace(σz() => -1, a'()*a() => 0)

# N-level constraint
P = Subspace(a'()*a() => 0)  # zero photons
```

### Classification Enums

The package uses classification enums for operator analysis:

```julia
DIAGONAL   # Operator preserves the subspace
RAISING    # Operator raises out of subspace P
LOWERING   # Operator lowers into subspace P
MIXED      # Operator has mixed character
```

---

## Internal Functions

These functions are not exported but may be useful for advanced users:

### Lie Algebra Support

```julia
# Get generators for detected Lie algebra
get_generators_for_lie_system(lie_info::NamedTuple)

# Convert between bases
gellmann_to_cartan_weyl(V_od, N, algebra_id)
cartan_weyl_to_gellmann(transitions, N, generators)
```

### Operator Classification

```julia
# Classify a single term
classify_term(term::QuTerm, coeff, constraints)

# Check if operator contains only specific types
has_only_bosons(term::QuTerm)
has_only_tls(term::QuTerm)
```
