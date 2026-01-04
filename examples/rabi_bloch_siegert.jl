#=
Rabi Model: Beyond Rotating Wave Approximation
===============================================

The quantum Rabi model describes a two-level system coupled to a harmonic 
oscillator without making the rotating wave approximation (RWA):

    H = ω a†a + Δ/2 σz + g σx(a + a†)
      = ω a†a + Δ/2 σz + g(σ⁺ + σ⁻)(a + a†)

This contains both:
- Jaynes-Cummings terms: g(a†σ⁻ + a σ⁺) - "rotating" terms
- Counter-rotating terms: g(a†σ⁺ + a σ⁻) - "anti-rotating" terms

The counter-rotating terms lead to the Bloch-Siegert shift, a correction
to the resonance frequency beyond the RWA.

In the dispersive regime with RWA, χ = g²/Δ
Including counter-rotating terms: χ_BS = g²/Δ + g²/(Δ+2ω)

Reference: 
- Bloch & Siegert, Phys. Rev. 57, 522 (1940)
- Forn-Díaz et al., Rev. Mod. Phys. 91, 025005 (2019)
=#

using UnitaryTransformations
using QuantumAlgebra
using Symbolics

println("="^60)
println("Rabi Model: Bloch-Siegert Shift")
println("="^60)

# Use σ± basis
QuantumAlgebra.use_σpm(true)

# Clear cached variables
UnitaryTransformations.clear_param_cache!()

# Define symbolic parameters
ω = Pr"ω"  # oscillator frequency
Δ = Pr"Δ"  # qubit splitting
g = Pr"g"  # coupling strength

# Full Rabi Hamiltonian (without RWA)
# H = ω a†a + Δ/2 σz + g(σ⁺ + σ⁻)(a + a†)
#   = ω a†a + Δ/2 σz + g(a†σ⁻ + a σ⁺) + g(a†σ⁺ + a σ⁻)
#                       ^^^^^^^^^^^^     ^^^^^^^^^^^^
#                       JC terms         Counter-rotating

H_full = ω * a'()*a() + Δ/2 * σz() + g * (σp() + σm()) * (a() + a'())

println("\n1. FULL RABI HAMILTONIAN")
println("-"^40)
println("H = ω a†a + Δ/2 σz + g(σ⁺ + σ⁻)(a + a†)")
println("\nExpanded:")
println("H = ", normal_form(H_full))

# For comparison, Jaynes-Cummings (RWA) Hamiltonian
H_JC = ω * a'()*a() + Δ/2 * σz() + g * (a'()*σm() + a()*σp())

println("\nJaynes-Cummings (RWA):")
println("H_JC = ", H_JC)

# Define subspace
P = Subspace(σz() => -1)

println("\n2. SUBSPACE DEFINITION")
println("-"^40)
println("P = ground state subspace: σz → -1")

# Decompose full Rabi Hamiltonian
H_d, H_od = decompose(H_full, P)

println("\n3. HAMILTONIAN DECOMPOSITION (Full Rabi)")
println("-"^40)
println("H_diagonal     = ", H_d)
println("H_off-diagonal = ", H_od)

# SW transformation for full Rabi model
println("\n4. SCHRIEFFER-WOLFF TRANSFORMATION (order 2)")
println("-"^40)
println("Transforming the full Rabi model...")

result_full = schrieffer_wolff(H_full, P; order=2)

println("Generator S = ", result_full.S)

# Also do JC for comparison
println("\nFor comparison, transforming JC model...")
result_JC = schrieffer_wolff(H_JC, P; order=2)

# Analyze effective Hamiltonians
println("\n5. EFFECTIVE HAMILTONIANS")
println("-"^40)

println("Full Rabi H_eff:")
terms_full = collect_terms(result_full.H_eff)
for (op, coeff) in terms_full
    println("  ", coeff, "  ", op)
end

println("\nJC (RWA) H_eff:")
terms_JC = collect_terms(result_JC.H_eff)
for (op, coeff) in terms_JC
    println("  ", coeff, "  ", op)
end

# Extract dispersive shifts
println("\n6. DISPERSIVE SHIFT COMPARISON")
println("-"^40)

χ_full = extract_coefficient(result_full.H_P, a'()*a())
χ_JC = extract_coefficient(result_JC.H_P, a'()*a())

println("Cavity frequency shift in ground state:")
println("  Full Rabi: χ = ", χ_full)
println("  JC (RWA):  χ = ", χ_JC)

println("\n7. BLOCH-SIEGERT SHIFT ANALYSIS")
println("-"^40)
println("""
The Bloch-Siegert shift is the correction from counter-rotating terms.

In RWA (Jaynes-Cummings):
  χ_JC = -g²/Δ   (for Δ = ω_q - ω_c > 0)

With counter-rotating terms, we expect additional contributions
from processes where both qubit and cavity are excited/de-excited together.

For near-resonance (Δ << ω), the counter-rotating terms contribute:
  χ_CR ≈ -g²/(Δ + 2ω)

The Bloch-Siegert shift is typically small when Δ << ω.
""")

# Numerical comparison
println("8. NUMERICAL VERIFICATION")
println("-"^40)
println("Parameters: ω = 5.0, Δ = 1.0, g = 0.1")
println("(Dispersive regime: g << Δ << ω)")

values = Dict(:ω => 5.0, :Δ => 1.0, :g => 0.1)

H_P_full_num = substitute_values(result_full.H_P, values)
H_P_JC_num = substitute_values(result_JC.H_P, values)

println("\nProjected Hamiltonians (numerical):")
println("  Full Rabi: H_P = ", H_P_full_num)
println("  JC (RWA):  H_P = ", H_P_JC_num)

# Calculate expected values
g_val = 0.1
Δ_val = 1.0
ω_val = 5.0

χ_JC_expected = -g_val^2 / Δ_val
χ_CR_expected = -g_val^2 / (Δ_val + 2*ω_val)  # Counter-rotating contribution
χ_full_expected = χ_JC_expected + χ_CR_expected

println("\nExpected dispersive shifts:")
println("  JC only: χ = -g²/Δ = ", χ_JC_expected)
println("  Counter-rotating: χ_CR = -g²/(Δ+2ω) = ", χ_CR_expected)
println("  Full (approx): χ_full ≈ ", χ_full_expected)
println("  Bloch-Siegert correction: ", χ_CR_expected, " (", 100*abs(χ_CR_expected/χ_JC_expected), "% of JC)")

println("\n9. PHYSICAL INTERPRETATION")
println("-"^40)
println("""
The counter-rotating terms (a†σ⁺ and a σ⁻) represent:
- Simultaneous excitation of qubit AND cavity (a†σ⁺)
- Simultaneous de-excitation of qubit AND cavity (a σ⁻)

These processes violate energy conservation by ~2ω and are suppressed,
but still contribute virtual corrections to the energy levels.

The Bloch-Siegert shift becomes important when:
1. The coupling g approaches the ultra-strong regime (g/ω > 0.1)
2. Precise frequency measurements are needed
3. Working with low-frequency oscillators

In circuit QED, this shift is typically ~1% or less of the JC dispersive shift.
""")
