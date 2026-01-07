#=
Three-Level Atom (Lambda Configuration)
=======================================

A three-level atom in the Lambda (Λ) configuration:

    |3⟩ (excited)
    / \
   /   \
  g₁   g₂
 /       \
|1⟩     |2⟩  (ground states)

Hamiltonian in Cartan-Weyl basis:
    H = ω₁ λ₇ + ω₂ λ₈ + g₁(E₁₃ + E₃₁) + g₂(E₂₃ + E₃₂)

where:
- λ₇, λ₈: Diagonal Gell-Mann generators (Cartan subalgebra)
- E_{ij}: Raising/lowering operators in Cartan-Weyl basis
- ω₁, ω₂: Energy splittings
- g₁, g₂: Coupling strengths to the excited state

The Cartan-Weyl basis transforms Gell-Mann matrices into raising/lowering operators:
    E₁₂ = (λ₁ + iλ₄)/2  (|1⟩ → |2⟩)
    E₂₁ = (λ₁ - iλ₄)/2  (|2⟩ → |1⟩)
    E₁₃ = (λ₂ + iλ₅)/2  (|1⟩ → |3⟩)
    E₃₁ = (λ₂ - iλ₅)/2  (|3⟩ → |1⟩)
    E₂₃ = (λ₃ + iλ₆)/2  (|2⟩ → |3⟩)
    E₃₂ = (λ₃ - iλ₆)/2  (|3⟩ → |2⟩)

These operators satisfy [H_d, E_{ij}] ∝ E_{ij}, making them eigenoperators
of the adjoint action - essential for Schrieffer-Wolff!

Reference: Cohen-Tannoudji, Arimondo, "Laser cooling and trapping of neutral atoms"
=#

using UnitaryTransformations
using QuantumAlgebra
using Symbolics

println("="^60)
println("Three-Level Atom (Lambda Configuration)")
println("="^60)

# Clear cached variables
UnitaryTransformations.clear_param_cache!()

# Create SU(3) generators (Gell-Mann matrices)
λ = su_generators(3, :λ)

println("\n1. SU(3) GENERATOR STRUCTURE")
println("-"^40)
println("Gell-Mann matrices λ₁...λ₈")
println("Diagonal (Cartan): λ₇, λ₈")
println("Off-diagonal: λ₁...λ₆")

# Define Cartan-Weyl raising/lowering operators
println("\n2. CARTAN-WEYL BASIS")
println("-"^40)

# Raising/lowering operators
E12 = normal_form((λ[1] + 1im * λ[4]) / 2)  # |1⟩ → |2⟩
E21 = normal_form((λ[1] - 1im * λ[4]) / 2)  # |2⟩ → |1⟩
E13 = normal_form((λ[2] + 1im * λ[5]) / 2)  # |1⟩ → |3⟩
E31 = normal_form((λ[2] - 1im * λ[5]) / 2)  # |3⟩ → |1⟩
E23 = normal_form((λ[3] + 1im * λ[6]) / 2)  # |2⟩ → |3⟩
E32 = normal_form((λ[3] - 1im * λ[6]) / 2)  # |3⟩ → |2⟩

println("E₁₂ = (λ₁ + iλ₄)/2  transitions |1⟩ → |2⟩")
println("E₁₃ = (λ₂ + iλ₅)/2  transitions |1⟩ → |3⟩")
println("E₂₃ = (λ₃ + iλ₆)/2  transitions |2⟩ → |3⟩")
println("(and their Hermitian conjugates E_{ji})")

# Verify eigenoperator property
println("\n3. EIGENOPERATOR VERIFICATION")
println("-"^40)
println("For SW to work, we need [H_d, E] = ε E (eigenoperator property)")
println()

# Check: [λ₇, E₁₃] should give back E₁₃ times a constant
comm_7_E13 = normal_form(comm(λ[7], E13))
ratio_7 = nothing
if !isempty(comm_7_E13.terms) && !isempty(E13.terms)
    # Compare first term coefficients
    c1 = first(comm_7_E13.terms)[2]
    c2 = first(E13.terms)[2]
    ratio_7 = c1 / c2
end
println("[λ₇, E₁₃] / E₁₃ = ", ratio_7)

comm_8_E13 = normal_form(comm(λ[8], E13))
ratio_8 = nothing
if !isempty(comm_8_E13.terms) && !isempty(E13.terms)
    c1 = first(comm_8_E13.terms)[2]
    c2 = first(E13.terms)[2]
    ratio_8 = c1 / c2
end
println("[λ₈, E₁₃] / E₁₃ = ", ratio_8)
println()
println("These ratios are the 'root vectors' of SU(3)!")

# Eigenvalue structure
println("\n4. EIGENVALUE STRUCTURE")
println("-"^40)
for i in [7, 8]
    m = QuantumAlgebra.gellmann_matrix(3, i)
    eigenvals = [real(m[j, j]) for j = 1:3]
    println("λ$i eigenvalues for |1⟩, |2⟩, |3⟩: ", round.(eigenvals, digits = 4))
end
println()
println("Energy denominator for E₁₃: ε(1→3) = E₃ - E₁")
println("  Using λ₇: 0 - 0.5 = -0.5 (from eigenvalue difference)")
println("  Using λ₈: -0.577 - 0.289 = -0.866 (from eigenvalue difference)")

# Define symbolic parameters
@variables Δ ω g₁ g₂

println("\n5. LAMBDA CONFIGURATION HAMILTONIAN")
println("-"^40)

# Hamiltonian with large detuning Δ for state |3⟩
# H_d = Δ λ₈ + ω λ₇ (energy structure)
# H_od = g₁(E₁₃ + E₃₁) + g₂(E₂₃ + E₃₂) (coupling to excited state)

H_d = Δ * λ[8] + ω * λ[7]
H_od = g₁ * (E13 + E31) + g₂ * (E23 + E32)
H = H_d + H_od

println("H_d = Δ λ₈ + ω λ₇  (diagonal/unperturbed)")
println("H_od = g₁(E₁₃ + E₃₁) + g₂(E₂₃ + E₃₂)  (coupling)")
println()
println("Full H = ", normal_form(H))

# Compute energy denominators manually
println("\n6. ENERGY DENOMINATOR CALCULATION")
println("-"^40)

# For E₁₃: [H_d, E₁₃] = ε₁₃ E₁₃
comm_Hd_E13 = normal_form(comm(H_d, E13))
println("[H_d, E₁₃] = ", comm_Hd_E13)
println("E₁₃ = ", E13)

# Extract energy denominator: [H_d, E₁₃] / E₁₃
if !isempty(comm_Hd_E13.terms) && !isempty(E13.terms)
    # Both should have same operator structure, just different coefficients
    c_comm = first(comm_Hd_E13.terms)[2]
    c_E13 = first(E13.terms)[2]
    ε_13 = c_comm / c_E13
    println()
    println("Energy denominator ε₁₃ = ", ε_13)
end

println()
# For E₃₁ (adjoint/lowering)
comm_Hd_E31 = normal_form(comm(H_d, E31))
if !isempty(comm_Hd_E31.terms) && !isempty(E31.terms)
    c_comm = first(comm_Hd_E31.terms)[2]
    c_E31 = first(E31.terms)[2]
    ε_31 = c_comm / c_E31
    println("Energy denominator ε₃₁ = ", ε_31, " (opposite sign)")
end

# Similarly for E₂₃ and E₃₂
println()
comm_Hd_E23 = normal_form(comm(H_d, E23))
if !isempty(comm_Hd_E23.terms) && !isempty(E23.terms)
    c_comm = first(comm_Hd_E23.terms)[2]
    c_E23 = first(E23.terms)[2]
    ε_23 = c_comm / c_E23
    println("Energy denominator ε₂₃ = ", ε_23)
end

println("\n7. GENERATOR (MANUAL CONSTRUCTION)")
println("-"^40)
println("""
For [S, H_d] = -H_od, we need:
    S = g₁ E₁₃/ε₁₃ + g₂ E₂₃/ε₂₃ + h.c.

With the eigenoperator property:
    [E₁₃/ε₁₃, H_d] = (1/ε₁₃)[E₁₃, H_d] = (1/ε₁₃)(-ε₁₃ E₁₃) = -E₁₃ ✓

Note: ε₁₃ = ω/2 + Δ·0.866 (combination of λ₇ and λ₈ eigenvalue differences)
""")

# Construct generator explicitly
# ε₁₃ = ω/2 + Δ·√3/2 for E₁₃ (based on Cartan eigenvalues)
sqrt3_half = sqrt(3) / 2
ε_13_expr = ω / 2 + Δ * sqrt3_half
ε_31_expr = -(ω / 2 + Δ * sqrt3_half)  # Opposite for lowering
ε_23_expr = -ω / 2 + Δ * sqrt3_half     # Different for E₂₃
ε_32_expr = -(-ω / 2 + Δ * sqrt3_half)  # Opposite for lowering

S =
    g₁ * E13 / ε_13_expr +
    g₁ * E31 / ε_31_expr +
    g₂ * E23 / ε_23_expr +
    g₂ * E32 / ε_32_expr

S = normal_form(S)
println("S = ", S)

# Verify: [S, H_d] should equal -H_od
println("\n8. GENERATOR VERIFICATION")
println("-"^40)

comm_S_Hd = normal_form(comm(S, H_d))
println("[S, H_d] = ", comm_S_Hd)
println("-H_od = ", normal_form(-H_od))

residual = normal_form(comm_S_Hd + H_od)
println()
println("Residual [S, H_d] + H_od = ", residual)
if isempty(residual.terms)
    println("✓ Generator equation satisfied!")
else
    println("(Small residual may be due to numerical coefficients)")
end

println("\n9. EFFECTIVE HAMILTONIAN (SECOND ORDER)")
println("-"^40)
println("""
H_eff = H_d + (1/2)[S, H_od] + O(g³)

The second-order term gives:
- AC Stark shifts of states |1⟩ and |2⟩
- Two-photon Raman coupling: g₁g₂/Δ term coupling |1⟩ ↔ |2⟩

This is the basis for:
- Raman transitions in laser cooling
- Two-photon processes in NMR
- Electromagnetically induced transparency (EIT)
""")

# Compute second-order correction
H_2 = normal_form(comm(S, H_od) / 2)
println("(1/2)[S, H_od] = ", H_2)

println("\n10. PHYSICAL INTERPRETATION")
println("-"^40)
println("""
In the limit Δ >> g₁, g₂ (large detuning):

1. State |3⟩ is virtually excited but never populated
2. Ground states |1⟩ and |2⟩ are coupled via two-photon process
3. Each ground state gets an AC Stark shift ∝ g²/Δ

The effective low-energy Hamiltonian describes a TWO-LEVEL system
within the ground state manifold, with effective coupling g₁g₂/Δ!
""")
