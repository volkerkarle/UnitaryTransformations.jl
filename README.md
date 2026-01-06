# UnitaryTransformations.jl

[![Build Status](https://github.com/volkerkarle/UnitaryTransformations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/volkerkarle/UnitaryTransformations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)

A Julia package for **symbolic unitary transformations** on quantum Hamiltonians. Built on [QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl) and [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/volkerkarle/UnitaryTransformations.jl")
```

## Quick Example

Derive the dispersive shift for qubit readout in circuit QED:

```julia
using UnitaryTransformations, QuantumAlgebra, Symbolics

QuantumAlgebra.use_σpm(true)
@variables Δ g

# Jaynes-Cummings Hamiltonian
H = Δ/2 * σz() + g * (a'()*σm() + a()*σp())

# Schrieffer-Wolff transformation
P = Subspace(σz() => -1)  # qubit ground state
result = schrieffer_wolff(H, P; order=2)

println(result.H_P)  # -Δ/2 + (-g²/Δ) a†a  ← dispersive shift χ = -g²/Δ
```

## Multi-Atom Systems with Symbolic Sums

For systems with many identical particles, use `SymSum` to represent symbolic sums that correctly generate **exchange interactions**:

```julia
using UnitaryTransformations, QuantumAlgebra, Symbolics
import QuantumAlgebra: sumindex, SymSum, SymExpr

QuantumAlgebra.use_σpm(true)
@variables ω_c Δ g

i = sumindex(1)

# Tavis-Cummings Hamiltonian: N atoms coupled to a cavity
H = SymExpr(ω_c * a'()*a()) + 
    SymSum(Δ/2 * σz(i), i) + 
    SymSum(g * (a'()*σm(i) + a()*σp(i)), i)

P = Subspace(a'()*a() => 0)  # Zero photon sector
result = schrieffer_wolff(H, P; order=2)

# The effective Hamiltonian includes exchange terms: χ Σᵢ≠ⱼ σ⁺ᵢσ⁻ⱼ
```

The `SymSum` type correctly handles:
- **Same-site terms**: `Σᵢ [Aᵢ, Bᵢ]`
- **Cross-site terms**: `Σᵢ Σⱼ≠ᵢ [Aᵢ, Bⱼ]` (exchange interactions!)

## Transformations

| Transformation | Purpose |
|----------------|---------|
| **Schrieffer-Wolff** | Block-diagonalize H, derive effective low-energy Hamiltonians |
| **Magnus Expansion** | Effective Hamiltonians for periodically driven (Floquet) systems |

## Supported Systems

| System | Operators |
|--------|-----------|
| Two-level systems | `σx()`, `σy()`, `σz()`, `σp()`, `σm()` |
| Bosonic modes | `a()`, `a'()` |
| N-level atoms | `nlevel_ops(N, :name)` |
| SU(N) algebras | `su_generators(N, :name)` |
| Multi-atom systems | `SymSum(expr, index)` for symbolic sums |
| Hybrid systems | Any combination |

## Features

- **Symbolic results**: Get analytical expressions like `-g²/Δ`, not floating-point numbers
- **Arbitrary perturbation order**: Compute to order 2, 4, 6+ with optional parallel acceleration
- **Automatic method selection**: Eigenoperator method for TLS/bosons, matrix-element method for SU(N)
- **Symbolic sums**: `SymSum` correctly handles multi-particle commutators with exchange terms

```julia
# Higher-order with parallelization
result = schrieffer_wolff(H, P; order=4, parallel=true)
```

## Magnus Expansion

For periodically driven systems H(t) = Σₙ Hₙ e^{inωt}:

```julia
@variables Δ Ω ω

modes = Dict(0 => Δ/2 * σz(), 1 => Ω/2 * σp(), -1 => Ω/2 * σm())
result = magnus_expansion(modes, ω; order=4)
# Computes Bloch-Siegert shift and higher-order corrections
```

## Documentation

- [Schrieffer-Wolff](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/schrieffer_wolff/) — Theory, tutorial, and examples
- [Magnus Expansion](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/magnus/) — Floquet systems and driven dynamics
- [API Reference](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/api/) — Function documentation

## Requirements

- Julia 1.12+
- [QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl)
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

## License

MIT

## Citation

```bibtex
@software{UnitaryTransformations.jl,
  author = {Karle, Volker},
  title = {UnitaryTransformations.jl: Symbolic Unitary Transformations for Quantum Hamiltonians},
  url = {https://github.com/volkerkarle/UnitaryTransformations.jl},
  year = {2025}
}
```

## Contributing

Contributions welcome! See the [documentation](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/) for details.
