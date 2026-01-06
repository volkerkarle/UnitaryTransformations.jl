# UnitaryTransformations.jl

A Julia package for performing symbolic unitary transformations on quantum Hamiltonians.

![Summary of Schrieffer-Wolff transformations at different orders](assets/summary.png)

## Overview

Unitary transformations are essential tools in quantum mechanics for simplifying Hamiltonians, deriving effective low-energy theories, and block-diagonalizing systems with separated energy scales.

This package provides **symbolic implementations** that produce analytical expressions rather than numerical results. For example, the dispersive shift in circuit QED is computed as `-g²/Δ`, not as a floating-point number.

## Available Transformations

| Transformation | Purpose | Documentation |
|----------------|---------|---------------|
| [**Schrieffer-Wolff**](@ref schrieffer_wolff) | Block-diagonalize Hamiltonians, derive effective low-energy theories | [Guide](@ref schrieffer_wolff) |
| [**Magnus Expansion**](@ref magnus_expansion) | Effective Hamiltonians for periodically driven (Floquet) systems | [Guide](@ref magnus_expansion) |

## Quick Start

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

QuantumAlgebra.use_σpm(true)

@variables Δ g

# Jaynes-Cummings Hamiltonian
H = Δ/2 * σz() + g * (a'()*σm() + a()*σp())

# Define the low-energy subspace
P = Subspace(σz() => -1)

# Perform Schrieffer-Wolff transformation
result = schrieffer_wolff(H, P; order=2)

println(result.H_P)
# Output: -Δ/2 + (-g²/Δ) a†a
```

## Supported Quantum Systems

The package works with quantum systems provided by [QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl):

| System | Operators | Example |
|--------|-----------|---------|
| **Two-level systems** | `σx()`, `σy()`, `σz()`, `σp()`, `σm()` | Qubits, spin-1/2 |
| **Bosonic modes** | `a()`, `a'()` | Cavities, phonons |
| **N-level atoms** | `nlevel_ops(N, :name)` | Multi-level atoms |
| **SU(N) systems** | `su_generators(N, :name)` | 3-level Λ systems |
| **Fermions** | `f(:name)`, `f'(:name)` | Electrons |

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/volkerkarle/UnitaryTransformations.jl")
```

## Documentation

```@contents
Pages = ["schrieffer_wolff.md", "magnus.md", "api.md"]
Depth = 2
```

## Citation

If you use this package in your research:

```bibtex
@software{UnitaryTransformations.jl,
  author = {Karle, Volker},
  title = {UnitaryTransformations.jl: Symbolic Unitary Transformations for Quantum Hamiltonians},
  url = {https://github.com/volkerkarle/UnitaryTransformations.jl},
  year = {2025}
}
```
