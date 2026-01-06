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

## Transformations

| Transformation | Purpose | Status |
|----------------|---------|--------|
| **Schrieffer-Wolff** | Block-diagonalize H, derive effective low-energy Hamiltonians | Implemented |
| **Magnus Expansion** | Effective Hamiltonians for periodically driven (Floquet) systems | Implemented |
| **Lang-Firsov** | Polaron transformation (electron-phonon coupling) | Planned |
| **Bogoliubov** | Diagonalize quadratic bosonic Hamiltonians | Planned |
| **Holstein-Primakoff** | Map spin operators to bosons | Planned |

## Supported Systems

| System | Operators |
|--------|-----------|
| Two-level systems | `σx()`, `σy()`, `σz()`, `σp()`, `σm()` |
| Bosonic modes | `a()`, `a'()` |
| N-level atoms | `nlevel_ops(N, :name)` |
| SU(N) algebras | `su_generators(N, :name)` |
| Hybrid systems | Any combination |

## Features

- **Symbolic results**: Get analytical expressions like `-g²/Δ`, not floating-point numbers
- **Arbitrary perturbation order**: Compute to order 2, 4, 6+ with optional parallel acceleration
- **Automatic method selection**: Eigenoperator method for TLS/bosons, matrix-element method for SU(N)

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

- [Tutorial](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/tutorial/) — Step-by-step introduction
- [Theory](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/theory/) — Mathematical foundations
- [Examples](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/examples/) — Physics applications
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
