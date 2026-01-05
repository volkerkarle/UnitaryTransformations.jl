# UnitaryTransformations.jl

A Julia package for performing symbolic unitary transformations on quantum Hamiltonians.

![Summary of Schrieffer-Wolff transformations at different orders](assets/summary.png)

## Overview

Unitary transformations are essential tools in quantum mechanics for:

- **Simplifying Hamiltonians** by eliminating unwanted couplings
- **Deriving effective theories** that capture low-energy physics
- **Block-diagonalizing** systems with separated energy scales
- **Changing to more convenient representations** (e.g., polaron frame)

This package provides symbolic implementations that produce analytical expressions rather than numerical results. For example, the dispersive shift in circuit QED is computed as `-g²/Δ`, not as a floating-point number.

## Available Transformations

| Transformation | Purpose | Status |
|----------------|---------|--------|
| [**Schrieffer-Wolff**](@ref sw_transformation) | Block-diagonalize Hamiltonians, derive effective low-energy theories | ✓ Implemented |
| **Lang-Firsov** | Eliminate linear electron-phonon coupling (polaron frame) | Planned |
| **Bogoliubov** | Diagonalize quadratic bosonic Hamiltonians | Planned |
| **Holstein-Primakoff** | Map spin operators to bosonic operators | Planned |

## Quick Start

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

# Use σ± basis (recommended)
QuantumAlgebra.use_σpm(true)

# Define symbolic parameters
@variables Δ g

# Jaynes-Cummings Hamiltonian
H = Δ/2 * σz() + g * (a'()*σm() + a()*σp())

# Define the low-energy subspace
P = Subspace(σz() => -1)  # qubit ground state

# Perform Schrieffer-Wolff transformation
result = schrieffer_wolff(H, P; order=2)

# The effective Hamiltonian
println(result.H_P)
# Output: -Δ/2 + (-g²/Δ) a†a
```

## Supported Quantum Systems

The package works with a variety of quantum systems provided by [QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl):

| System | Operators | Example |
|--------|-----------|---------|
| **Two-level systems** | `σx()`, `σy()`, `σz()`, `σp()`, `σm()` | Qubits, spin-1/2 |
| **Bosonic modes** | `a()`, `a'()` | Cavities, phonons |
| **N-level atoms** | `nlevel_ops(N, :name)` | Multi-level atoms |
| **SU(N) systems** | `su_generators(N, :name)` | 3-level Λ systems |
| **Fermions** | `f(:name)`, `f'(:name)` | Electrons |

These can be combined freely—for example, an N-level atom coupled to a bosonic cavity.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/volkerkarle/UnitaryTransformations.jl")
```

Dependencies (automatically installed):
- [QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl) — Symbolic quantum operator algebra
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) — Computer algebra system

## Documentation Structure

```@contents
Pages = ["theory.md", "tutorial.md", "examples.md", "api.md"]
Depth = 1
```

### [Theory](@ref theory)
Mathematical background on unitary transformations, including the Schrieffer-Wolff method and the Baker-Campbell-Hausdorff formula.

### [Tutorial](@ref tutorial)
Step-by-step guide to using the Schrieffer-Wolff transformation.

### [Examples](@ref examples)
Physics applications: Jaynes-Cummings, Rabi model, multi-level atoms.

### [API Reference](@ref api)
Complete function documentation.

## Package Philosophy

1. **Symbolic over numeric**: Results are analytical expressions that can be simplified, manipulated, and substituted.

2. **Automatic method selection**: The package chooses the optimal algorithm based on the operator types present.

3. **Extensible design**: New transformations can be added as separate modules while sharing the infrastructure.

4. **Verified correctness**: Extensive tests verify mathematical identities and compare with known results.

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
