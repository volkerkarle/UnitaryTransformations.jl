# [Examples](@id examples)

This page presents complete physics examples demonstrating the Schrieffer-Wolff transformation. Each example shows the LaTeX output for easy use in publications.

## Convergence Analysis

The following figure demonstrates how higher-order SW transformations systematically improve the accuracy of the effective Hamiltonian:

![Summary of SW transformations at different orders](assets/summary.png)

**Key observations:**
- **Panel (a)**: Energy levels of a two-level system. Higher orders capture the curvature better.
- **Panel (b)**: The energy gap ``E_+ - E_-`` approaches the exact ``\sqrt{\Delta^2 + 4\varepsilon^2}`` with higher orders.
- **Panel (c)**: Relative error decreases with increasing SW order.
- **Panel (d)**: Jaynes-Cummings model showing dispersive shift ``\chi`` and Kerr nonlinearity ``K`` (only present at order 4+).

## Jaynes-Cummings: Dispersive Regime

The Jaynes-Cummings model describes a two-level atom coupled to a single cavity mode:

```math
H = \omega_c a^\dagger a + \frac{\omega_q}{2}\sigma_z + g(a^\dagger\sigma^- + a\sigma^+)
```

In the **dispersive regime** (``|\Delta| = |\omega_q - \omega_c| \gg g``), the Schrieffer-Wolff transformation yields an effective Hamiltonian with a state-dependent frequency shift.

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

QuantumAlgebra.use_σpm(true)

@variables Δ g  # Δ = detuning, g = coupling

# Hamiltonian (in rotating frame)
H = Δ/2 * σz() + g * (a'()*σm() + a()*σp())

# Ground state subspace
P = Subspace(σz() => -1)

# Transform and display results
result = schrieffer_wolff(H, P; order=2)
show_result(result)
```

**Generator:**
```math
S = \frac{-g}{2\Delta} a^{\dagger} \sigma^{-} + \frac{g}{2\Delta} \sigma^{+} a
```

**Effective Hamiltonian:**
```math
H_{\text{eff}} = -\frac{\Delta}{2} - \frac{g^{2}}{\Delta} a^{\dagger} a + \frac{g^{2} + \Delta^{2}}{\Delta} \sigma^{+} \sigma^{-} + \frac{2g^{2}}{\Delta} a^{\dagger} \sigma^{+} \sigma^{-} a
```

**Projected to subspace P** (qubit ground state):
```math
H_P = -\frac{\Delta}{2} - \frac{g^{2}}{\Delta} a^{\dagger} a
```

### Physical Interpretation

The effective Hamiltonian contains:

1. **Dispersive shift**: ``\chi = -g^2/\Delta`` — cavity frequency shifts when qubit is in ground state
2. **AC Stark shift**: Qubit frequency shifts with photon number

This is the basis for **qubit readout** in circuit QED!

![Jaynes-Cummings dispersive shift and Kerr nonlinearity](assets/jaynes_cummings.png)

### Extracting Parameters

```julia
χ = extract_coefficient(result.H_P, a'()*a())
println(to_latex(χ))  # Output: \frac{-g^{2}}{\Delta}

# Numerical evaluation
H_num = substitute_values(result.H_P, Dict(:g => 0.1, :Δ => 1.0))
# χ = -0.01, matching -g²/Δ
```

---

## Two-Level System with Transverse Field

A qubit in longitudinal and transverse magnetic fields:

```math
H = \frac{\Delta}{2}\sigma_z + \varepsilon\sigma_x
```

This textbook problem has an exact solution, making it perfect for verification.

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

QuantumAlgebra.use_σpm(true)

@variables Δ ε

H = Δ/2 * σz() + ε * (σp() + σm())
P = Subspace(σz() => -1)

result = schrieffer_wolff(H, P; order=2)
print_latex(result.H_P; name="H_P")
```

**Output:**
```math
H_P = -\frac{\Delta}{2} - \frac{\varepsilon^{2}}{\Delta}
```

### Comparison with Exact Solution

The exact ground state energy is:

```math
E_- = -\sqrt{\frac{\Delta^2}{4} + \varepsilon^2} \approx -\frac{\Delta}{2} - \frac{\varepsilon^2}{\Delta} + O(\varepsilon^4)
```

The SW result matches the perturbation expansion exactly!

![Two-level system: SW orders vs exact solution](assets/two_level_system.png)

| ``\varepsilon/\Delta`` | Exact | SW (2nd order) | Error |
|:----------------------:|:-----:|:--------------:|:-----:|
| 0.01 | -0.50005 | -0.5001 | 0.00% |
| 0.05 | -0.50125 | -0.5025 | 0.00% |
| 0.10 | -0.50499 | -0.51 | 0.02% |
| 0.20 | -0.51980 | -0.54 | 0.28% |

### Convergence with Order

![Convergence of SW expansion](assets/convergence.png)

The figure shows how the approximation error decreases with increasing SW order. For small perturbations (``\varepsilon/\Delta < 0.3``), convergence is rapid. For larger coupling, more orders are needed.

---

## Rabi Model: Bloch-Siegert Shift

The **full Rabi model** includes counter-rotating terms neglected in the rotating-wave approximation:

```math
H = \omega a^\dagger a + \frac{\Delta}{2}\sigma_z + g(\sigma^+ + \sigma^-)(a + a^\dagger)
```

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

QuantumAlgebra.use_σpm(true)

@variables ω Δ g

H = ω * a'()*a() + Δ/2 * σz() + g * (σp() + σm()) * (a() + a'())
P = Subspace(σz() => -1)

result = schrieffer_wolff(H, P; order=2)
print_latex(result.H_P; name="H_P")
```

**Output:**
```math
H_P = -\frac{\Delta}{2} + \omega a^{\dagger} a - \frac{g^{2}}{\Delta - \omega} a^{\dagger} a - \frac{g^{2}}{\Delta + \omega} a^{\dagger} a + \frac{g^{2}}{\Delta + \omega} - \frac{g^{2}}{\Delta + \omega} (a^{\dagger})^2 - \frac{g^{2}}{\Delta + \omega} a^2
```

### Physical Interpretation

Compared to Jaynes-Cummings, the full Rabi model produces:

1. **JC dispersive shift**: ``-g^2/(\Delta - \omega)`` from rotating terms
2. **Bloch-Siegert shift**: ``-g^2/(\Delta + \omega)`` from counter-rotating terms
3. **Squeezing terms**: ``a^2`` and ``(a^\dagger)^2`` that squeeze the cavity field

The total dispersive shift combines both contributions:

```math
\chi_{\text{total}} = -\frac{g^2}{\Delta - \omega} - \frac{g^2}{\Delta + \omega}
```

![Bloch-Siegert shift: JC vs full Rabi model](assets/bloch_siegert.png)

---

## N-Level Atom in a Cavity

For atoms with more than two levels, use `nlevel_ops`:

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

# 5-level atom: σ[i,j] = |i⟩⟨j|
σ5 = nlevel_ops(5, :q)

ω = [Symbolics.variable(Symbol("ω", i)) for i in 1:5]
@variables ωc g

# Atom + cavity + dipole coupling |1⟩↔|3⟩
H = sum(ω[i] * σ5[i,i] for i in 1:5) + 
    ωc * a'()*a() + 
    g * (σ5[1,3] * a'() + σ5[3,1] * a())

# Zero-photon subspace
P = Subspace(a'()*a() => 0)

result = schrieffer_wolff(H, P; order=2)
show_result(result)
```

**Generator:**
```math
S = \frac{g}{\omega_1 - \omega_3 - \omega_c} |1\rangle\langle 3| \, a^{\dagger} + \frac{g}{\omega_3 - \omega_1 + \omega_c} |3\rangle\langle 1| \, a
```

### Physical Interpretation

- **Dispersive shift**: ``\chi_{13} = g^2/(\omega_1 - \omega_3 + \omega_c)``
- **AC Stark shifts** on levels 1 and 3
- **Other levels** (2, 4, 5) appear only with bare energies

---

## Three-Level Lambda System

For systems with SU(N) symmetry, use `su_generators`:

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

# Gell-Mann matrices for SU(3)
λ = su_generators(3, :λ)

@variables Δ Ω₁ Ω₂

# Lambda system
H = Δ * λ[8] + Ω₁ * λ[1] + Ω₂ * λ[4]

P = Subspace(λ[8] => 1/sqrt(3))

result = schrieffer_wolff(H, P; order=2)
show_result(result)
```

The package automatically detects SU(3) and uses the matrix-element method.

### When to Use SU(N) vs N-Level

| Approach | Use When |
|----------|----------|
| `nlevel_ops` | Physical atoms, specific transitions |
| `su_generators` | Systems with SU(N) symmetry |

---

## Running the Examples

Complete example files are in the `examples/` directory:

```bash
julia --project examples/jaynes_cummings_dispersive.jl
julia --project examples/two_level_system.jl
julia --project examples/rabi_bloch_siegert.jl
julia --project examples/three_level_atom.jl
```

---

## Tips

1. **Use `show_result(result)`** to see all components in LaTeX
2. **Use `to_latex(expr)`** to get a LaTeX string for any expression
3. **Use `extract_coefficient(expr, op)`** to get specific parameters
4. **Use `substitute_values(expr, Dict(...))`** for numerical evaluation
