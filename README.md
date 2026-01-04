# UnitaryTransformations.jl

[![Build Status](https://github.com/volkerkarle/UnitaryTransformations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/volkerkarle/UnitaryTransformations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)
[![Julia 1.12+](https://img.shields.io/badge/Julia-1.12+-blue.svg)](https://julialang.org/)

A Julia package for performing symbolic Schrieffer-Wolff transformations on quantum Hamiltonians.

## What is this?

**Schrieffer-Wolff transformation** is a perturbative method to block-diagonalize quantum Hamiltonians with well-separated energy scales. Given a Hamiltonian `H = H₀ + V` where `V` couples different energy sectors, the transformation finds an effective Hamiltonian `H_eff` that acts only within a chosen subspace.

This package provides:
- **Symbolic computation** of effective Hamiltonians with proper energy denominators
- **Automatic handling** of commutator algebra for bosons, fermions, and spins
- **Integration with Symbolics.jl** for clean expressions like `g²/Δ`

Built on [QuantumAlgebra.jl](https://github.com/volkerkarle/QuantumAlgebra.jl) for symbolic quantum operator algebra.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/volkerkarle/UnitaryTransformations.jl")
```

## Quick Start

### Jaynes-Cummings Dispersive Shift

The classic application: a qubit coupled to a cavity in the dispersive regime.

```julia
using UnitaryTransformations, QuantumAlgebra, Symbolics

# Use σ± basis for clean results
QuantumAlgebra.use_σpm(true)

# Define symbolic parameters
@variables Δ g  # Δ = qubit-cavity detuning, g = coupling strength

# Define Hamiltonian: H = Δ/2 σz + g(a†σ⁻ + a σ⁺)
H = Δ/2 * σz() + g * (a'()*σm() + a()*σp())

# Define subspace: qubit in ground state
P = Subspace(σz() => -1)

# Perform Schrieffer-Wolff to second order
result = schrieffer_wolff(H, P; order=2)

# Effective Hamiltonian projected to ground state subspace
println(result.H_P)
# Output: -0.5Δ + (-g²/Δ) a†a
#                  ^^^^^^
#         This is the dispersive shift χ = -g²/Δ
```

### Two-Level System

A spin in longitudinal and transverse fields:

```julia
using UnitaryTransformations, QuantumAlgebra, Symbolics

QuantumAlgebra.use_σpm(true)

# Define symbolic parameters
@variables Δ ε

# H = Δ/2 σz + ε σx  (where σx = σ⁺ + σ⁻)
H = Δ/2 * σz() + ε * (σp() + σm())
P = Subspace(σz() => -1)

result = schrieffer_wolff(H, P; order=2)

# Ground state energy: E = -Δ/2 - ε²/Δ
# (matches 2nd-order perturbation theory)
println(collect_terms(result.H_P))
```

## Key Features

### Symbolic Energy Denominators

Unlike numerical approaches, this package produces symbolic results with proper energy denominators:

```julia
# The dispersive shift is computed symbolically
χ = extract_coefficient(result.H_P, a'()*a())
# Returns: -(g^2) / Δ
```

### Simplification and Substitution

```julia
# Simplify coefficients
H_simplified = simplify_coefficients(result.H_eff)

# Substitute numerical values
H_numeric = substitute_values(result.H_P, Dict(:g => 0.1, :Δ => 1.0))
```

### Verified Against Exact Solutions

The package has been rigorously tested:
- Generator equation `[S, H_d] = -V_od` verified to be exactly satisfied
- Numerical results match exact eigenvalues to <0.02% for ε/Δ = 0.1
- Dispersive shift correctly gives `-g²/Δ`

## Examples

See the `examples/` folder for complete worked examples:

- **`jaynes_cummings_dispersive.jl`** - Qubit-cavity dispersive shift (circuit QED)
- **`two_level_system.jl`** - Spin in transverse field with exact solution comparison  
- **`rabi_bloch_siegert.jl`** - Full Rabi model including counter-rotating terms

## API Overview

### Core Functions

| Function | Description |
|----------|-------------|
| `schrieffer_wolff(H, P; order=2)` | Main SW transformation |
| `Subspace(σz() => -1)` | Define low-energy subspace |
| `decompose(H, P)` | Split H into diagonal/off-diagonal parts |
| `solve_for_generator(H_d, V_od, P)` | Find generator S |

### Result Structure

`schrieffer_wolff` returns a named tuple with:
- `H_eff` - Block-diagonal effective Hamiltonian
- `S` - Generator of the transformation
- `H_P` - H_eff projected onto subspace P

### Utility Functions

| Function | Description |
|----------|-------------|
| `simplify_coefficients(expr)` | Simplify Symbolics coefficients |
| `substitute_values(expr, Dict(...))` | Substitute numerical values |
| `extract_coefficient(expr, op)` | Extract coefficient of operator |
| `collect_terms(expr)` | List all terms with coefficients |

## Theory

The Schrieffer-Wolff transformation finds a unitary `U = e^S` such that:

```
H_eff = e^S H e^{-S} = H + [S,H] + ½[S,[S,H]] + ...
```

is block-diagonal with respect to a chosen subspace P. The generator S is determined by solving:

```
[S, H_d] = -V_od
```

where `H_d` is the diagonal part and `V_od` is the off-diagonal perturbation.

## Requirements

- Julia 1.12+
- QuantumAlgebra.jl (automatically installed)
- Symbolics.jl (automatically installed)

## License

GPL-3.0

## Citation

If you use this package in your research, please cite:

```bibtex
@software{UnitaryTransformations.jl,
  author = {Karle, Volker},
  title = {UnitaryTransformations.jl: Symbolic Schrieffer-Wolff Transformations},
  url = {https://github.com/volkerkarle/UnitaryTransformations.jl},
  year = {2024}
}
```
