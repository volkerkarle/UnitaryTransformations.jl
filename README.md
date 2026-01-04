# UnitaryTransformations.jl

[![Build Status](https://github.com/volkerkarle/UnitaryTransformations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/volkerkarle/UnitaryTransformations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)

A Julia package for performing **symbolic unitary transformations** on quantum Hamiltonians.

## Overview

Unitary transformations are fundamental tools in quantum mechanics for simplifying Hamiltonians, eliminating unwanted couplings, and deriving effective low-energy theories. This package provides symbolic implementations of several important transformations:

| Transformation | Purpose | Status |
|----------------|---------|--------|
| **Schrieffer-Wolff** | Block-diagonalize H, derive effective low-energy Hamiltonians | Implemented |
| **Lang-Firsov** | Eliminate linear electron-phonon coupling (polaron transformation) | Planned |
| **Bogoliubov** | Diagonalize quadratic bosonic Hamiltonians | Planned |
| **Holstein-Primakoff** | Map spin operators to bosonic operators | Planned |

Built on [QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl) for symbolic quantum operator algebra and [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) for symbolic mathematics.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/volkerkarle/UnitaryTransformations.jl")
```

## Quick Example: Dispersive Readout

The Schrieffer-Wolff transformation derives the dispersive Hamiltonian used for qubit readout in circuit QED:

```julia
using UnitaryTransformations, QuantumAlgebra, Symbolics

QuantumAlgebra.use_σpm(true)
@variables Δ g  # Δ = detuning, g = coupling

# Jaynes-Cummings Hamiltonian
H = Δ/2 * σz() + g * (a'()*σm() + a()*σp())

# Transform to eliminate qubit-cavity coupling
P = Subspace(σz() => -1)  # qubit ground state
result = schrieffer_wolff(H, P; order=2)

println(result.H_P)
# Output: -Δ/2 + (-g²/Δ) a†a
#                 ^^^^^^ dispersive shift χ = -g²/Δ
```

## Supported Quantum Systems

The package handles a wide range of quantum systems:

| System Type | Operators | Example Use |
|-------------|-----------|-------------|
| **Two-level systems** | `σx()`, `σy()`, `σz()`, `σp()`, `σm()` | Qubits, spin-1/2 |
| **Bosonic modes** | `a()`, `a'()` | Cavities, phonons |
| **N-level atoms** | `nlevel_ops(N, :name)` | Multi-level atoms, transmons |
| **SU(N) algebras** | `su_generators(N, :name)` | Collective spin, Lambda systems |
| **Hybrid systems** | Combinations | Atom + cavity, spin + phonon |

### Example: Multi-Level Atom + Cavity

```julia
using UnitaryTransformations, QuantumAlgebra, Symbolics

# 5-level atom with transition operators σ[i,j] = |i⟩⟨j|
σ5 = nlevel_ops(5, :q)

# Symbolic level energies
ω = [Symbolics.variable(Symbol("ω", i)) for i in 1:5]
@variables ωc g

# Atom + cavity + dipole coupling between levels 1↔3
H = sum(ω[i] * σ5[i,i] for i in 1:5) + ωc * a'()*a() + 
    g * (σ5[1,3] * a'() + σ5[3,1] * a())

# SW transformation in cavity vacuum
P = Subspace(a'()*a() => 0)
result = schrieffer_wolff(H, P; order=2)

# Result: dispersive shifts, AC Stark corrections
```

## Key Features

### Symbolic Results

Unlike numerical approaches, this package produces **analytical expressions**:

```julia
χ = extract_coefficient(result.H_P, a'()*a())
# Returns: -g²/Δ  (symbolic, not floating-point!)
```

### Automatic Method Selection

The package automatically chooses the optimal algorithm based on operator types:

| Operator Type | Method | Description |
|---------------|--------|-------------|
| TLS, bosons, N-level | Eigenoperator | Uses `[H₀, O] = εO` structure |
| SU(N) generators | Matrix-element | Works in Cartan-Weyl basis |

### Arbitrary Perturbation Order

```julia
result_2nd = schrieffer_wolff(H, P; order=2)
result_4th = schrieffer_wolff(H, P; order=4)  # Higher accuracy
```

### Verified Implementation

- Generator equation `[S, H₀] = -V` verified exactly
- Numerical accuracy <0.02% in perturbative regime
- 147 tests covering TLS, bosons, SU(3), and N-level systems

## Documentation

- **[Tutorial](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/tutorial/)**: Step-by-step introduction
- **[Theory](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/theory/)**: Mathematical foundations (BCH, SW derivation)
- **[Examples](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/examples/)**: Physics applications (JC, Rabi, N-level)
- **[API Reference](https://volkerkarle.github.io/UnitaryTransformations.jl/dev/api/)**: Function documentation

## Transformations

### Schrieffer-Wolff (Implemented)

Finds a unitary U = eˢ that block-diagonalizes H = H₀ + V:

```
H_eff = eˢ H e⁻ˢ = H₀ + ½[S, V] + O(V³)
```

The generator S satisfies `[S, H₀] = -V`. 

**Applications**:
- Dispersive readout in circuit QED
- Exchange interactions (t-J model from Hubbard)
- Effective spin models
- Adiabatic elimination of fast modes

### Lang-Firsov (Planned)

Eliminates linear electron-phonon coupling via phonon displacement:

```
U = exp(∑ₖ gₖ/ωₖ (bₖ† - bₖ) n)
```

**Applications**: Polarons, phonon sidebands, Franck-Condon physics

### Bogoliubov (Planned)

Diagonalizes quadratic bosonic Hamiltonians:

```
H = ∑ₖ ωₖ aₖ†aₖ + Δₖ(aₖa₋ₖ + h.c.)  →  H = ∑ₖ Ωₖ βₖ†βₖ
```

**Applications**: BCS superconductivity, squeezed states, magnons

### Holstein-Primakoff (Planned)

Maps spin-S operators to bosons for large-S expansions:

```
S⁺ → √(2S) √(1 - a†a/2S) a
```

**Applications**: Spin waves, magnon-polaritons

## Requirements

- Julia 1.10+
- [QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl) (with SU(N) support)
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

Contributions are welcome! Particularly for:
- New transformation types (Lang-Firsov, Bogoliubov, Holstein-Primakoff)
- Additional physics examples
- Performance improvements
- Documentation enhancements

Please open an issue to discuss before submitting a PR.
