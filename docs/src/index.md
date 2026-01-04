# UnitaryTransformations.jl

A Julia package for performing symbolic Schrieffer-Wolff transformations on quantum Hamiltonians.

## Overview

The **Schrieffer-Wolff transformation** is a powerful perturbative technique for deriving effective low-energy Hamiltonians. Given a quantum system with well-separated energy scales, it systematically eliminates high-energy degrees of freedom while capturing their effects through renormalized parameters.

This package provides:

- **Symbolic computation** of effective Hamiltonians
- **Automatic commutator algebra** for bosons, fermions, and spins  
- **Proper energy denominators** like `g²/Δ` using Symbolics.jl
- **Projection** onto chosen low-energy subspaces

### When to Use Schrieffer-Wolff

The transformation is ideal when:

1. Your Hamiltonian has a **small perturbation** coupling different energy sectors
2. You want to work in a **reduced subspace** (e.g., ground state manifold)
3. You need **analytical expressions** for effective parameters

Common applications include:
- Dispersive readout in circuit QED
- Exchange interactions from virtual hopping (t-J model from Hubbard)
- Effective spin models from electronic systems

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/volkerkarle/UnitaryTransformations.jl")
```

The package will automatically install its dependencies:
- [QuantumAlgebra.jl](https://github.com/volkerkarle/QuantumAlgebra.jl) - Symbolic quantum operator algebra
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) - Computer algebra system

## Quick Start

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

# Use σ± basis (recommended for SW transformations)
QuantumAlgebra.use_σpm(true)

# Define symbolic parameters
@variables Δ g  # Δ = energy splitting, g = coupling strength

# Jaynes-Cummings Hamiltonian
H = Δ/2 * σz() + g * (a'()*σm() + a()*σp())

# Define the low-energy subspace (qubit in ground state)
P = Subspace(σz() => -1)

# Perform Schrieffer-Wolff transformation
result = schrieffer_wolff(H, P; order=2)

# The effective Hamiltonian in the ground state subspace
println(result.H_P)
# Output: -0.5Δ + (-(g^2))/Δ a†a
```

## Package Structure

```
UnitaryTransformations.jl/
├── src/
│   ├── UnitaryTransformations.jl  # Main module
│   ├── subspace.jl                # Subspace type definition
│   ├── decompose.jl               # H → H_d + V_od decomposition
│   ├── commutator_series.jl       # BCH expansion
│   ├── inverse_liouvillian.jl     # Solve [S, H_d] = -V_od
│   ├── schrieffer_wolff.jl        # Main SW algorithm
│   └── symbolic_utils.jl          # Simplification utilities
├── examples/
│   ├── jaynes_cummings_dispersive.jl
│   ├── two_level_system.jl
│   └── rabi_bloch_siegert.jl
└── test/
    └── runtests.jl
```

## Next Steps

- See the [Tutorial](@ref tutorial) for a step-by-step guide
- Check the [Examples](@ref examples) for physics applications
- Browse the [API Reference](@ref api) for function documentation
