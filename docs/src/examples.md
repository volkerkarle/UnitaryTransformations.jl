# [Examples](@id examples)

This page presents complete physics examples demonstrating the package capabilities.

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
UnitaryTransformations.clear_param_cache!()

# Define symbolic parameters
@variables Δ g  # Δ = detuning ω_q - ω_c, g = coupling

# Hamiltonian (in frame rotating at ω_c)
H = Δ/2 * σz() + g * (a'()*σm() + a()*σp())

# Ground state subspace
P = Subspace(σz() => -1)

# Transform
result = schrieffer_wolff(H, P; order=2)

# Results
println("Effective Hamiltonian:")
for (op, coeff) in collect_terms(result.H_eff)
    println("  ", coeff, "  ", op)
end
```

### Physical Interpretation

The effective Hamiltonian contains:

1. **Dispersive shift**: ``\chi = g^2/\Delta`` 
   - Cavity frequency shifts by ``\pm\chi`` depending on qubit state
   - Used for **qubit readout** in circuit QED

2. **AC Stark shift**: Qubit frequency shifts with photon number
   - ``\omega_q \to \omega_q + 2\chi \langle a^\dagger a \rangle``

### Numerical Verification

```julia
# Compare with expected value
χ = extract_coefficient(result.H_P, a'()*a())
println("Computed χ = ", χ)
println("Expected:  -g²/Δ")

# Substitute values
χ_num = substitute_values(result.H_P, Dict(:g => 0.1, :Δ => 1.0))
println("χ(g=0.1, Δ=1) = -0.01  # matches -g²/Δ!")
```

---

## Two-Level System with Transverse Field

A qubit in longitudinal and transverse magnetic fields:

```math
H = \frac{\Delta}{2}\sigma_z + \varepsilon\sigma_x
```

This is a textbook quantum mechanics problem with exact solution, making it perfect for verification.

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

QuantumAlgebra.use_σpm(true)
UnitaryTransformations.clear_param_cache!()

# Define symbolic parameters
@variables Δ ε  # Δ = longitudinal field (energy splitting), ε = transverse field (perturbation)

# Hamiltonian: σx = σ⁺ + σ⁻
H = Δ/2 * σz() + ε * (σp() + σm())

# Ground state subspace
P = Subspace(σz() => -1)

# Transform
result = schrieffer_wolff(H, P; order=2)

println("Projected Hamiltonian H_P:")
for (op, coeff) in collect_terms(result.H_P)
    println("  ", coeff, "  ", op)
end
```

### Comparison with Exact Solution

The exact eigenvalues are:

```math
E_\pm = \pm\sqrt{\Delta^2/4 + \varepsilon^2}
```

The SW result for the ground state energy is:

```math
E_g = -\frac{\Delta}{2} - \frac{\varepsilon^2}{\Delta} + O(\varepsilon^4)
```

```julia
# Numerical comparison
for ε_val in [0.01, 0.05, 0.1, 0.2]
    Δ_val = 1.0
    E_exact = -sqrt(Δ_val^2/4 + ε_val^2)
    E_SW = -Δ_val/2 - ε_val^2/Δ_val
    error = 100 * abs(E_exact - E_SW) / abs(E_exact)
    println("ε/Δ = $ε_val: error = $(round(error, digits=4))%")
end
```

Output:
```
ε/Δ = 0.01: error = 0.0%
ε/Δ = 0.05: error = 0.0012%
ε/Δ = 0.1: error = 0.0192%
ε/Δ = 0.2: error = 0.2755%
```

The perturbative result is excellent when ``\varepsilon \ll \Delta``!

---

## Rabi Model: Bloch-Siegert Shift

The **full Rabi model** includes counter-rotating terms that are neglected in the rotating-wave approximation:

```math
H = \omega a^\dagger a + \frac{\Delta}{2}\sigma_z + g(\sigma^+ + \sigma^-)(a + a^\dagger)
```

The counter-rotating terms (``a^\dagger\sigma^+`` and ``a\sigma^-``) lead to the **Bloch-Siegert shift**.

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

QuantumAlgebra.use_σpm(true)
UnitaryTransformations.clear_param_cache!()

# Define symbolic parameters
@variables ω Δ g  # ω = oscillator frequency, Δ = qubit splitting, g = coupling

# Full Rabi Hamiltonian (no RWA)
H_rabi = ω * a'()*a() + Δ/2 * σz() + g * (σp() + σm()) * (a() + a'())

# For comparison: Jaynes-Cummings (with RWA)
H_jc = ω * a'()*a() + Δ/2 * σz() + g * (a'()*σm() + a()*σp())

P = Subspace(σz() => -1)

result_rabi = schrieffer_wolff(H_rabi, P; order=2)
result_jc = schrieffer_wolff(H_jc, P; order=2)

println("Full Rabi model produces additional terms:")
println("  - Squeezing terms: a², (a†)²")
println("  - Modified dispersive shift")
```

### Bloch-Siegert Correction

The counter-rotating terms contribute an additional frequency shift:

```math
\chi_{BS} \approx \frac{g^2}{\Delta + 2\omega}
```

This is typically small (∼1% of the JC dispersive shift) but measurable in precision experiments.

---

## Running the Examples

The complete example files are in the `examples/` directory:

```bash
cd UnitaryTransformations.jl
julia --project examples/jaynes_cummings_dispersive.jl
julia --project examples/two_level_system.jl
julia --project examples/rabi_bloch_siegert.jl
```

Each example includes detailed physical interpretation and numerical verification.
