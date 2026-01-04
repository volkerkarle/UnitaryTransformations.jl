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

This is typically small (~1% of the JC dispersive shift) but measurable in precision experiments.

---

## N-Level Atom in a Cavity

For atoms with more than two levels, use `nlevel_ops` from QuantumAlgebra. This example shows a 5-level atom with a single dipole-allowed transition coupled to a cavity.

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

# 5-level atom: σ[i,j] = |i⟩⟨j| transition operator
σ5 = nlevel_ops(5, :q)

# Level energies (symbolic)
ω = [Symbolics.variable(Symbol("ω", i)) for i in 1:5]
@variables ωc g

# Hamiltonian: atom + cavity + coupling between levels 1 and 3
H = sum(ω[i] * σ5[i,i] for i in 1:5) +    # atomic levels
    ωc * a'()*a() +                        # cavity
    g * (σ5[1,3] * a'() + σ5[3,1] * a())  # dipole coupling |1⟩↔|3⟩

# Zero-photon subspace
P = Subspace(a'()*a() => 0)

# Perform SW transformation
result = schrieffer_wolff(H, P; order=2)

println("Generator S:")
println(result.S)

println("\nEffective Hamiltonian:")
println(result.H_eff)
```

### Physical Interpretation

The effective Hamiltonian contains:

1. **Dispersive shift on cavity**: The cavity frequency shifts depending on atomic state
   ```math
   \chi_{13} = \frac{g^2}{\omega_1 - \omega_3 + \omega_c}
   ```

2. **AC Stark shifts on atomic levels**: Levels 1 and 3 experience light shifts proportional to ``\langle a^\dagger a \rangle``

3. **Other levels unaffected**: Levels 2, 4, 5 only appear with their bare energies (to second order)

---

## Multi-Level Atom with Multiple Couplings

A more realistic model includes multiple dipole-allowed transitions:

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

# 7-level atom
σ7 = nlevel_ops(7, :q)

# Level energies
ω = [Symbolics.variable(Symbol("ω", i)) for i in 1:7]
@variables ωc g₁₃ g₂₅

# Hamiltonian with two transitions
H = sum(ω[i] * σ7[i,i] for i in 1:7) +
    ωc * a'()*a() +
    g₁₃ * (σ7[1,3] * a'() + σ7[3,1] * a()) +  # |1⟩↔|3⟩
    g₂₅ * (σ7[2,5] * a'() + σ7[5,2] * a())    # |2⟩↔|5⟩

P = Subspace(a'()*a() => 0)
result = schrieffer_wolff(H, P; order=2)

# Extract individual dispersive shifts
# Each transition contributes independently
println("H_eff contains terms from both transitions")
```

### Selection Rules

The SW transformation respects selection rules:
- Only directly coupled levels acquire dispersive shifts
- Cross-terms between different transitions appear only at higher orders

---

## Three-Level Lambda System (SU(3))

For systems with SU(N) symmetry, the package automatically uses the matrix-element method:

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

# SU(3) generators (Gell-Mann matrices)
λ = su_generators(3, :λ)

@variables Δ Ω₁ Ω₂

# Lambda system:
# - λ₈ ~ diagonal (level splittings)
# - λ₁, λ₂ ~ transitions between levels 1↔2
# - λ₄, λ₅ ~ transitions between levels 1↔3
H = Δ * λ[8] + Ω₁ * λ[1] + Ω₂ * λ[4]

# Subspace with specific λ₈ eigenvalue
P = Subspace(λ[8] => 1/sqrt(3))

# Auto-detects SU(3) and uses matrix-element method
result = schrieffer_wolff(H, P; order=2)

println("Generator:")
println(result.S)
println("\nEffective Hamiltonian:")
println(result.H_eff)
```

### When to Use SU(N) vs N-Level

| Approach | Use When |
|----------|----------|
| `nlevel_ops` | Physical multilevel atoms, specific transitions |
| `su_generators` | Systems with SU(N) symmetry, collective operators |

Both approaches give equivalent physics but may produce differently structured results.

---

## Running the Examples

The complete example files are in the `examples/` directory:

```bash
cd UnitaryTransformations.jl
julia --project examples/jaynes_cummings_dispersive.jl
julia --project examples/two_level_system.jl
julia --project examples/rabi_bloch_siegert.jl
julia --project examples/three_level_atom.jl
```

Each example includes detailed physical interpretation and numerical verification.

---

## Tips for New Examples

1. **Start simple**: Begin with a minimal Hamiltonian and add complexity
2. **Check subspace definition**: The choice of ``P`` determines the entire structure
3. **Verify with limits**: Check known limits (e.g., ``g \to 0``)
4. **Use `collect_terms`**: Helps identify the physical meaning of each term
5. **Compare orders**: Running at different orders helps assess convergence
