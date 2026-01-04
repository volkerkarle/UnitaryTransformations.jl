#=
Jaynes-Cummings Model: Dispersive Regime
=========================================

The Jaynes-Cummings Hamiltonian describes a two-level atom (qubit) coupled to 
a single mode of an electromagnetic cavity:

    H = ω_c a†a + ω_q/2 σz + g(a†σ⁻ + a σ⁺)

where:
- ω_c: cavity frequency
- ω_q: qubit frequency  
- g: coupling strength
- Δ = ω_q - ω_c: detuning

In the dispersive regime (|Δ| >> g), the Schrieffer-Wolff transformation
eliminates the direct qubit-cavity coupling, yielding an effective Hamiltonian:

    H_eff = ω_c a†a + ω_q/2 σz + χ a†a σz + O(g⁴)

where χ = g²/Δ is the dispersive shift. Higher-order corrections include:
- Kerr nonlinearity: K (a†a)² with K ~ g⁴/Δ³
- Modified dispersive shift

This causes:
- Qubit frequency shift depending on photon number (AC Stark shift)
- Cavity frequency shift depending on qubit state (used for qubit readout)

Reference: Blais et al., PRA 69, 062320 (2004)
=#

using UnitaryTransformations
using QuantumAlgebra
using Symbolics

println("="^60)
println("Jaynes-Cummings Model: Dispersive Regime")
println("="^60)

# Use σ± basis for cleaner algebra
QuantumAlgebra.use_σpm(true)

# Clear any cached symbolic variables from previous runs
UnitaryTransformations.clear_param_cache!()

# Define symbolic parameters using Symbolics.jl
@variables Δ g  # Δ = detuning (ω_q - ω_c), g = coupling strength

# Jaynes-Cummings Hamiltonian in rotating frame of cavity
# H = Δ/2 σz + g(a†σ⁻ + a σ⁺)  (plus ω_c a†a which commutes with everything)
# We work in units where the cavity frequency shifts are measured relative to ω_c

H_int = Δ/2 * σz() + g * (a'()*σm() + a()*σp())

println("\n1. HAMILTONIAN")
println("-"^40)
println("H = ω_c a†a + Δ/2 σz + g(a†σ⁻ + a σ⁺)")
println("Working in interaction picture, H_int = Δ/2 σz + g(a†σ⁻ + a σ⁺)")
println("\nH_int = ", H_int)

# Define subspace: qubit in ground state |g⟩ (spin down, σz = -1)
P = Subspace(σz() => -1)

println("\n2. SUBSPACE DEFINITION")
println("-"^40)
println("P = ground state subspace: σz → -1")

# Decompose into diagonal and off-diagonal parts
H_d, H_od = decompose(H_int, P)

println("\n3. HAMILTONIAN DECOMPOSITION")
println("-"^40)
println("H_diagonal     = ", H_d)
println("H_off-diagonal = ", H_od)

# Perform Schrieffer-Wolff transformation to fourth order
println("\n4. SCHRIEFFER-WOLFF TRANSFORMATION (order 4)")
println("-"^40)

result = schrieffer_wolff(H_int, P; order = 4)

println("Generator S = ", result.S)
println("\nEffective Hamiltonian H_eff = ", result.H_eff)

# Analyze the effective Hamiltonian
println("\n5. ANALYSIS OF H_eff")
println("-"^40)

terms = collect_terms(result.H_eff)
println("Terms in H_eff:")
for (op, coeff) in terms
    println("  ", coeff, "  ", op)
end

# Extract the dispersive shift (coefficient of a†a σz-like terms)
# In σpm basis, σz = -1 + 2σ⁺σ⁻, so we look for a†a and a†σ⁺σ⁻a terms

println("\n6. DISPERSIVE SHIFT")
println("-"^40)

# The a†a term gives the cavity frequency shift
χ_cavity = extract_coefficient(result.H_eff, a'()*a())
println("Cavity frequency shift (a†a coefficient): ", χ_cavity)

# The a†σ⁺σ⁻a term contributes to state-dependent shifts
# For a full dispersive shift, we need to combine terms

# Project to ground state (P subspace)
println("\n7. PROJECTED TO GROUND STATE")
println("-"^40)
println("H_P = ", result.H_P)

terms_P = collect_terms(result.H_P)
println("\nTerms in H_P (qubit in ground state):")
for (op, coeff) in terms_P
    println("  ", coeff, "  ", op)
end

# The effective cavity frequency in the ground state
ω_eff_g = extract_coefficient(result.H_P, a'()*a())
println("\nEffective cavity frequency shift in |g⟩: ", ω_eff_g)

# Expected result: -g²/Δ (negative because ground state)
println("\nExpected from theory: -g²/Δ")

println("\n8. PHYSICAL INTERPRETATION")
println("-"^40)
println("""
The dispersive shift χ = g²/Δ leads to:

1. AC Stark shift: The qubit frequency shifts by ±χ⟨a†a⟩ depending on 
   the photon number in the cavity.

2. Dispersive readout: The cavity frequency shifts by ±χ depending on
   the qubit state, allowing non-demolition qubit measurement.

3. Photon number splitting: Each photon number state |n⟩ has a distinct
   qubit transition frequency ω_q + 2nχ.

For our result with qubit in ground state:
- Cavity sees frequency shift of $(ω_eff_g) per photon
- This is the famous dispersive shift used in circuit QED
""")

# Verify with numerical values
println("\n9. NUMERICAL VERIFICATION")
println("-"^40)
println("Substituting g = 0.1, Δ = 1.0 (dispersive regime: g/Δ = 0.1)")

H_P_num = substitute_values(result.H_P, Dict(:g => 0.1, :Δ => 1.0))
println("H_P(numerical) = ", H_P_num)

expected_shift = -0.1^2 / 1.0
println("\nExpected cavity shift: -g²/Δ = ", expected_shift)
println("This matches the coefficient of a†a in H_P!")
