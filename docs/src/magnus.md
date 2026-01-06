# [Magnus Expansion](@id magnus_expansion)

The **Magnus expansion** is a technique for solving time-dependent Schrödinger equations and computing effective Hamiltonians for periodically driven (Floquet) systems.

---

## Theory

### The Problem

Consider a time-dependent Schrödinger equation:

```math
i\frac{\partial U}{\partial t} = H(t) U(t)
```

where ``H(t)`` is periodic with period ``T = 2\pi/\omega``. We want to find an effective time-independent Hamiltonian ``H_{\text{eff}}`` such that after one period:

```math
U(T) = e^{-i H_{\text{eff}} T}
```

### Fourier Representation

A periodic Hamiltonian can be written in Fourier form:

```math
H(t) = \sum_n H_n e^{in\omega t}
```

where the Hermiticity condition requires ``H_{-n} = H_n^\dagger``.

### The Magnus Series

The solution to the Schrödinger equation can be written as:

```math
U(t) = e^{\Omega(t)}
```

where ``\Omega(t)`` is the Magnus series:

```math
\Omega(t) = \Omega_1(t) + \Omega_2(t) + \Omega_3(t) + \cdots
```

The effective Hamiltonian is:

```math
H_{\text{eff}} = i\Omega(T)/T = \sum_{k=1}^\infty \Omega_k
```

### Orders of the Expansion

**First order (k=1):**
```math
\Omega_1 = H_0
```

The leading term is simply the time-averaged Hamiltonian.

**Second order (k=2):**
```math
\Omega_2 = \sum_{n>0} \frac{-[H_n, H_{-n}]}{n\omega}
```

This captures effects like the **Bloch-Siegert shift** from counter-rotating drive terms.

**Higher orders (k≥3):**

For order ``k``, the Magnus term involves nested commutators with ``k`` Fourier components:

```math
\Omega_k = \sum_{\substack{n_1,\ldots,n_k \\ \sum_j n_j = 0}} C(n_1,\ldots,n_k) \, [[\cdots[[H_{n_1}, H_{n_2}], H_{n_3}],\ldots], H_{n_k}]
```

The coefficients are:

```math
C(n_1,\ldots,n_k) = \frac{1}{\omega^{k-1} \prod_{j=1}^{k-1} s_j}
```

where ``s_j = n_1 + \cdots + n_j`` are partial sums.

### Resonance Condition

The sum ``\sum_j n_j = 0`` is the **resonance condition**. It ensures that only certain combinations of Fourier modes contribute to the effective Hamiltonian.

### Reducible Terms

A term is **reducible** if any intermediate partial sum ``s_j = 0``. Such terms are excluded because they factorize into products of lower-order terms.

### Convergence

The Magnus expansion converges when:

```math
\int_0^T \|H(t)\| \, dt < \pi
```

For high-frequency driving (``\omega \gg \|H\|``), convergence is typically rapid.

### Physical Applications

| Application | Key Effect |
|-------------|------------|
| NMR | Rotating-wave corrections |
| Trapped ions | Micromotion effects |
| Circuit QED | AC Stark shifts |
| Floquet engineering | Artificial gauge fields |

### Comparison with Floquet Theory

The Magnus expansion provides a **high-frequency expansion** of Floquet theory. While exact Floquet diagonalization gives the full quasienergy spectrum, the Magnus expansion produces analytical expressions valid in the ``\omega \gg \|V\|`` regime.

---

## How-To Guide

This section walks through using the Magnus expansion step by step.

### Step 1: Define Fourier Modes

The Magnus expansion works with Hamiltonians in Fourier representation:

```math
H(t) = \sum_n H_n e^{in\omega t}
```

Provide the Fourier modes as a dictionary mapping mode index to operator:

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

QuantumAlgebra.use_σpm(true)

@variables Δ Ω ω

# H(t) = (Δ/2)σz + (Ω/2)(e^{iωt}σ⁺ + e^{-iωt}σ⁻)
modes = Dict(
    0  => Δ/2 * σz(),      # Time-independent part
    1  => Ω/2 * σp(),      # e^{iωt} coefficient
    -1 => Ω/2 * σm()       # e^{-iωt} coefficient
)
```

### Step 2: Compute the Expansion

Call `magnus_expansion` with the Fourier modes and driving frequency:

```julia
result = magnus_expansion(modes, ω; order=4)
```

This returns a named tuple with:
- `result.H_eff` - The total effective Hamiltonian
- `result.Ω1`, `result.Ω2`, ... - Individual order contributions
- `result.orders` - Dictionary of all computed orders

### Step 3: Analyze Results

```julia
println("H_eff = ", result.H_eff)
println("Ω₁ (time-averaged) = ", result.Ω1)
println("Ω₂ (Bloch-Siegert) = ", result.Ω2)
println("Ω₃ = ", result.Ω3)
```

### Hermiticity Check

The Magnus expansion verifies that your Hamiltonian is Hermitian:

```julia
# This will check H₋ₙ = Hₙ† automatically
result = magnus_expansion(modes, ω; check_hermitian=true)

# Skip the check if you know what you're doing
result = magnus_expansion(modes, ω; check_hermitian=false)
```

---

## Examples

### Circularly Driven Two-Level System

A qubit driven by a circularly polarized field:

```math
H(t) = \frac{\Delta}{2}\sigma_z + \frac{\Omega}{2}\left(e^{i\omega t}\sigma^+ + e^{-i\omega t}\sigma^-\right)
```

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

QuantumAlgebra.use_σpm(true)

@variables Δ Ω ω

# Fourier modes: H(t) = Σₙ Hₙ exp(inωt)
modes = Dict(
    0  => Δ/2 * σz(),
    1  => Ω/2 * σp(),
    -1 => Ω/2 * σm()
)

# Compute Magnus expansion to 4th order
result = magnus_expansion(modes, ω; order=4)

println("H_eff = ", result.H_eff)
println("Ω₂ (Bloch-Siegert) = ", result.Ω2)
println("Ω₃ = ", result.Ω3)
```

**Output (order 2):**
```math
\Omega_2 = -\frac{\Omega^2}{4\omega}\sigma_z
```

This is the **Bloch-Siegert shift** — a correction to the qubit frequency due to counter-rotating drive terms.

### Higher-Order Corrections

The Magnus expansion can be computed to arbitrary order:

```julia
# High-order expansion
result = magnus_expansion(modes, ω; order=6)

# Access individual orders
println("Order 3: ", result.Ω3)  # ~ -ΔΩ²/(4ω²)
println("Order 4: ", result.Ω4)  # ~ +Δ²Ω²/(4ω³)
println("Order 5: ", result.Ω5)  # ~ -Δ³Ω²/(4ω⁴)
```

**Pattern:** For k ≥ 2, each order contributes:

```math
\Omega_k = (-1)^{k-1} \frac{\Delta^{k-2} \Omega^2}{4\omega^{k-1}} \sigma_z
```

This geometric series sums to the exact Bloch-Siegert result:

```math
H_{\text{eff}} = \frac{\Delta}{2}\left(1 - \frac{\Omega^2}{2\omega(\omega + \Delta)}\right)\sigma_z
```

### Multiple Driving Frequencies

For Hamiltonians with multiple Fourier components:

```julia
@variables Δ Ω₁ Ω₂ ω

# Bichromatic driving
modes = Dict(
    0  => Δ/2 * σz(),
    1  => Ω₁/2 * σp(),
    -1 => Ω₁/2 * σm(),
    2  => Ω₂/2 * σp(),
    -2 => Ω₂/2 * σm()
)

result = magnus_expansion(modes, ω; order=3)
```

### Example: AC Stark Shift

When a qubit is driven off-resonantly, the Magnus expansion captures the AC Stark shift:

```julia
@variables δ Ω ω  # δ = detuning from qubit frequency

modes = Dict(
    0  => δ/2 * σz(),
    1  => Ω/2 * σp(),
    -1 => Ω/2 * σm()
)

result = magnus_expansion(modes, ω; order=2)

# The AC Stark shift appears in Ω₂
println("AC Stark shift: ", result.Ω2)
```

---

## Running the Examples

Complete example files are in the `examples/` directory:

```bash
julia --project examples/driven_qubit.jl  # Magnus expansion
```

---

## Key Points

1. **Provide Fourier modes** as a `Dict{Int, Operator}` mapping ``n \to H_n``
2. **Hermiticity** requires ``H_{-n} = H_n^\dagger`` — the package checks this by default
3. **Access individual orders** via `result.Ω2`, `result.Ω3`, etc.
4. **High-frequency limit**: The expansion is most accurate when ``\omega \gg \|H\|``
5. **Use `order` parameter** to control the number of terms computed

---

## References

1. W. Magnus, "On the exponential solution of differential equations for a linear operator," *Comm. Pure Appl. Math.* **7**, 649 (1954).

2. S. Blanes et al., "The Magnus expansion and some of its applications," *Physics Reports* **470**, 151 (2009).

3. A. Eckardt and E. Anisimovas, "High-frequency approximation for periodically driven quantum systems from a Floquet-space perspective," *New J. Phys.* **17**, 093039 (2015).

4. M. Bukov, L. D'Alessio, and A. Polkovnikov, "Universal high-frequency behavior of periodically driven systems: from dynamical stabilization to Floquet engineering," *Adv. Phys.* **64**, 139 (2015).
