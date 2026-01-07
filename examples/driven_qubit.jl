# Driven Qubit Example: Magnus Expansion
# 
# This example demonstrates the Magnus expansion for a driven two-level system,
# showing how to compute effective Hamiltonians for periodic driving.

using QuantumAlgebra
using UnitaryTransformations
using Symbolics

# Use σ± basis for cleaner results
QuantumAlgebra.use_σpm(true)

# Define symbolic parameters
@variables Δ Ω ω

println("="^60)
println("Magnus Expansion for Driven Two-Level Systems")
println("="^60)

# =============================================================================
# Example 1: Linearly polarized drive (σx coupling)
# =============================================================================

println("\n" * "="^60)
println("Example 1: Linearly polarized drive")
println("H(t) = Δ/2 σz + Ω cos(ωt) σx")
println("="^60)

# cos(ωt) = (e^{iωt} + e^{-iωt})/2
# σx = σ⁺ + σ⁻
# So H₁ = H₋₁ = Ω/2 × σx/2 = Ω/4 (σ⁺ + σ⁻) ... wait that's wrong
# Actually: Ω cos(ωt) σx = Ω/2 (e^{iωt} + e^{-iωt}) σx
# So H₁ = Ω/2 σx and H₋₁ = Ω/2 σx

modes_linear =
    Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * (σp() + σm()), -1 => Ω / 2 * (σp() + σm()))

println("\nFourier modes:")
println("  H₀ = Δ/2 σz")
println("  H₁ = H₋₁ = Ω/2 σx")

result_linear = magnus_expansion(modes_linear, ω; order = 2)

println("\nMagnus expansion results (order 2):")
println("  Ω₁ = ", result_linear.Ω1)
println("  Ω₂ = ", result_linear.Ω2)
println("  H_eff = ", result_linear.H_eff)

println("\nNote: Ω₂ = 0 because [σx, σx] = 0 (operator commutes with itself)")

# =============================================================================
# Example 2: Circularly polarized drive (Bloch-Siegert shift)
# =============================================================================

println("\n" * "="^60)
println("Example 2: Circularly polarized drive")
println("H(t) = Δ/2 σz + Ω(e^{iωt} σ⁺ + e^{-iωt} σ⁻)")
println("="^60)

modes_circular = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())

println("\nFourier modes:")
println("  H₀ = Δ/2 σz")
println("  H₁ = Ω/2 σ⁺")
println("  H₋₁ = Ω/2 σ⁻")

println("\nComputing Magnus expansion to order 3...")
result_circ = magnus_expansion(modes_circular, ω; order = 3)

println("\nMagnus expansion results:")
println("  Ω₁ = ", result_circ.Ω1)
println("  Ω₂ = ", result_circ.Ω2)
println("  Ω₃ = ", result_circ.Ω3)
println("\n  H_eff (total) = ", result_circ.H_eff)

println("\nPhysical interpretation:")
println("  - Ω₂ ~ Ω²/ω: AC Stark shift (Bloch-Siegert-like)")
println("  - Ω₃ ~ ΔΩ²/ω²: Higher-order correction")

# =============================================================================
# Example 3: Bichromatic drive (two frequencies)
# =============================================================================

println("\n" * "="^60)
println("Example 3: Bichromatic drive")
println("H(t) = Δ/2 σz + Ω₁ cos(ωt) σx + Ω₂ cos(2ωt) σx")
println("="^60)

@variables Ω₁ Ω₂

σx_expr = σp() + σm()

modes_bichromatic = Dict(
    0 => Δ / 2 * σz(),
    1 => Ω₁ / 2 * σx_expr,
    -1 => Ω₁ / 2 * σx_expr,
    2 => Ω₂ / 2 * σx_expr,
    -2 => Ω₂ / 2 * σx_expr,
)

println("\nFourier modes:")
println("  H₀ = Δ/2 σz")
println("  H₁ = H₋₁ = Ω₁/2 σx (fundamental)")
println("  H₂ = H₋₂ = Ω₂/2 σx (second harmonic)")

result_bichromatic = magnus_expansion(modes_bichromatic, ω; order = 2)

println("\nMagnus expansion results (order 2):")
println("  Ω₂ = ", result_bichromatic.Ω2)
println("  H_eff = ", result_bichromatic.H_eff)

# =============================================================================
# Example 4: Comparing with rotating frame approximation
# =============================================================================

println("\n" * "="^60)
println("Example 4: Comparison with RWA (Rotating Wave Approximation)")
println("="^60)

println("\nIn a frame rotating at frequency ω, the RWA gives:")
println("  H_RWA = (Δ-ω)/2 σz + Ω/2 σx")
println("")
println("The Magnus expansion captures corrections to RWA when ω ≉ Δ,")
println("including counter-rotating terms that RWA neglects.")
println("")
println("For Δ << ω (far-detuned drive), Magnus gives effective shifts")
println("proportional to Ω²/ω, similar to the Bloch-Siegert shift.")

# =============================================================================
# Numerical example
# =============================================================================

println("\n" * "="^60)
println("Numerical Example: Ω/ω = 0.1, Δ/ω = 0.5")
println("="^60)

# Substitute numerical values
using Symbolics: substitute

values = Dict(Δ => 0.5, Ω => 0.1, ω => 1.0)

H_eff_numeric = substitute(result_circ.H_eff, values)
println("\nH_eff (numeric) = ", H_eff_numeric)

# Also show individual contributions
Ω2_numeric = substitute(result_circ.Ω2, values)
Ω3_numeric = substitute(result_circ.Ω3, values)
println("Ω₂ contribution = ", Ω2_numeric)
println("Ω₃ contribution = ", Ω3_numeric)

println("\nRelative size of corrections:")
println("  |Ω₂/Ω₁| ~ (Ω/ω)² = ", 0.1^2)
println("  |Ω₃/Ω₁| ~ Δ(Ω/ω)²/ω = ", 0.5 * 0.1^2)

QuantumAlgebra.use_σpm(false)

println("\n" * "="^60)
println("Done!")
println("="^60)
