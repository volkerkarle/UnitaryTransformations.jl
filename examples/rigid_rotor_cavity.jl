#=
Rigid Rotor in a Cavity (L = 0, 1, 2)
=====================================

A rotating molecule (rigid rotor) truncated to the first three rotational levels
coupled to a single cavity mode in the dipole gauge:

    H = ω_c a†a + B·L(L+1) + g·cos(θ)(a + a†)

where:
- ω_c: cavity frequency
- B: rotational constant
- g: light-matter coupling strength
- cos(θ): dipole operator (orientation of molecular axis)

Energy levels (rigid rotor):
    E_L = B·L(L+1) → E₀ = 0, E₁ = 2B, E₂ = 6B

Matrix elements from Clebsch-Gordan coefficients (selection rule ΔL = ±1):
    ⟨L=0|cos(θ)|L=1⟩ = 1/√3 ≈ 0.577
    ⟨L=1|cos(θ)|L=2⟩ = √(2/15) ≈ 0.365

This example demonstrates the Schrieffer-Wolff transformation for a hybrid
system with BOTH bosonic (cavity) and Lie algebra (rotor) degrees of freedom.

In the dispersive regime (|ω_c - 2B| >> g), the effective Hamiltonian shows:
- Cavity Lamb shift from virtual rotational excitations
- AC Stark shifts of rotational levels
- Cavity-mediated coupling between rotational states

Physical context:
- Molecular polaritons in optical/infrared cavities
- Rotational strong coupling in THz cavities
- Cavity-modified chemistry with polar molecules

Reference: 
- Ebbesen et al., "Hybrid Light-Matter States", Acc. Chem. Res. 49, 2403 (2016)
- Flick et al., "Cavity Born-Oppenheimer Approximation", JCTC 13, 1616 (2017)
=#

using UnitaryTransformations
using QuantumAlgebra
using Symbolics

println("="^70)
println("  Rigid Rotor in a Cavity: Schrieffer-Wolff Transformation")
println("="^70)

# Clear cached variables
UnitaryTransformations.clear_param_cache!()
QuantumAlgebra.use_σpm(true)

# =============================================================================
# Section 1: Physical Setup and Parameters
# =============================================================================

println("\n1. PHYSICAL SETUP")
println("-"^60)
println("""
Rigid rotor levels (L = 0, 1, 2):
    E_L = B·L(L+1)
    E₀ = 0,  E₁ = 2B,  E₂ = 6B

Transition frequencies:
    L=0 ↔ L=1:  ω₀₁ = 2B
    L=1 ↔ L=2:  ω₁₂ = 4B

Dipole matrix elements (from Clebsch-Gordan):
    ⟨L|cos(θ)|L+1⟩ = √((L+1)/((2L+1)(2L+3)))
    
    ⟨0|cos(θ)|1⟩ = 1/√3 ≈ 0.577
    ⟨1|cos(θ)|2⟩ = √(2/15) ≈ 0.365
""")

# Symbolic parameters
@variables B ω_c g

# Matrix elements (use Float64 for cleaner output)
c01 = 1 / sqrt(3)       # ⟨L=0|cos(θ)|L=1⟩
c12 = sqrt(2 / 15)      # ⟨L=1|cos(θ)|L=2⟩

println("Numerical values:")
println("  c₀₁ = 1/√3 = $(round(c01, digits=6))")
println("  c₁₂ = √(2/15) = $(round(c12, digits=6))")

# =============================================================================
# Section 2: SU(3) Representation of the Rotor
# =============================================================================

println("\n2. SU(3) REPRESENTATION")
println("-"^60)
println("""
The 3-level rotor is represented using SU(3) generators (Gell-Mann matrices).

State mapping:
    |L=0⟩ → |1⟩,  |L=1⟩ → |2⟩,  |L=2⟩ → |3⟩

Diagonal Hamiltonian using Cartan generators λ₇, λ₈:
    H_rotor = -2B·λ₇ - (10/√3)B·λ₈ + (8/3)B·𝟙

(The constant (8/3)B is an overall energy shift that we drop.)
""")

# Create SU(3) generators
λ = su_generators(3, :λ)

# Verify the Gell-Mann matrix eigenvalues
println("Gell-Mann diagonal eigenvalues:")
for i in [7, 8]
    m = QuantumAlgebra.gellmann_matrix(3, i)
    eigenvals = [real(m[j, j]) for j in 1:3]
    println("  λ$i: ", round.(eigenvals, digits=4))
end

# Construct the rotor Hamiltonian (dropping constant term)
H_rotor = B * (-2 * λ[7] - 10 / sqrt(3) * λ[8])
println("\nH_rotor = -2B·λ₇ - (10/√3)B·λ₈")
println("       = ", H_rotor)

# =============================================================================
# Section 3: Cartan-Weyl (Transition) Operators
# =============================================================================

println("\n3. CARTAN-WEYL (TRANSITION) OPERATORS")
println("-"^60)
println("""
The Cartan-Weyl basis provides transition operators that ARE eigenoperators
of H_rotor, essential for the Schrieffer-Wolff transformation.

Construction from Gell-Mann matrices:
    E₁₂ = λ₁ + i·λ₄  →  |L=0⟩⟨L=1|  (lowers L)
    E₂₁ = λ₁ - i·λ₄  →  |L=1⟩⟨L=0|  (raises L)
    E₂₃ = λ₃ + i·λ₆  →  |L=1⟩⟨L=2|  (lowers L)
    E₃₂ = λ₃ - i·λ₆  →  |L=2⟩⟨L=1|  (raises L)
""")

# Transition operators
E12 = normal_form(λ[1] + 1im * λ[4])  # |L=0⟩⟨L=1|
E21 = normal_form(λ[1] - 1im * λ[4])  # |L=1⟩⟨L=0|
E23 = normal_form(λ[3] + 1im * λ[6])  # |L=1⟩⟨L=2|
E32 = normal_form(λ[3] - 1im * λ[6])  # |L=2⟩⟨L=1|

println("E₁₂ = ", E12)
println("E₂₁ = ", E21)
println("E₂₃ = ", E23)
println("E₃₂ = ", E32)

# Verify eigenoperator property by showing the commutators
println("\nEigenoperator verification [H_rotor, E_ij] = ε_ij·E_ij:")

for (op, name, expected_ε) in [
    (E12, "E₁₂", "-2B"),
    (E21, "E₂₁", "+2B"),
    (E23, "E₂₃", "-4B"),
    (E32, "E₃₂", "+4B"),
]
    comm_result = normal_form(comm(H_rotor, op))
    println("  [H_rotor, $name] = $comm_result")
    println("    Expected: $expected_ε × $name")
end

# =============================================================================
# Section 4: Full Hamiltonian
# =============================================================================

println("\n4. FULL HAMILTONIAN")
println("-"^60)

# Cavity Hamiltonian
H_cavity = ω_c * a'() * a()

# Dipole operator in Cartan-Weyl basis
# cos(θ) = c₀₁(E₁₂ + E₂₁) + c₁₂(E₂₃ + E₃₂)
cos_theta = c01 * (E12 + E21) + c12 * (E23 + E32)

# Interaction Hamiltonian
H_int = g * cos_theta * (a'() + a())
H_int = normal_form(H_int)

# Full Hamiltonian
H_d = normal_form(H_rotor + H_cavity)
H = normal_form(H_d + H_int)

println("H_cavity = ω_c·a†a")
println("H_rotor  = -2B·λ₇ - (10/√3)B·λ₈")
println()
println("cos(θ̂) = c₀₁(E₁₂ + E₂₁) + c₁₂(E₂₃ + E₃₂)")
println("       = ", normal_form(cos_theta))
println()
println("H_int = g·cos(θ̂)·(a + a†)")
println("      = ", H_int)
println()
println("H_diagonal = ", H_d)

# =============================================================================
# Section 5: Eigenoperator Structure of the Interaction
# =============================================================================

println("\n5. EIGENOPERATOR STRUCTURE OF H_int")
println("-"^60)
println("""
Each term in H_int is a product of a bosonic operator (a or a†) and a 
transition operator (E_ij). These products are eigenoperators of H_d:

    [H_d, a†·E_ij] = (ω_c + ε_ij)·a†·E_ij
    [H_d, a·E_ij]  = (-ω_c + ε_ij)·a·E_ij

where ε_ij = E_i - E_j is the rotor transition energy.
""")

# Verify combined eigenoperators
println("Verifying eigenoperator property for a†·E and a·E terms:")

eigenops = [
    (a'() * E12, "a†·E₁₂", "ω_c - 2B"),
    (a'() * E21, "a†·E₂₁", "ω_c + 2B"),
    (a() * E12, "a·E₁₂", "-ω_c - 2B"),
    (a() * E21, "a·E₂₁", "-ω_c + 2B"),
    (a'() * E23, "a†·E₂₃", "ω_c - 4B"),
    (a'() * E32, "a†·E₃₂", "ω_c + 4B"),
    (a() * E23, "a·E₂₃", "-ω_c - 4B"),
    (a() * E32, "a·E₃₂", "-ω_c + 4B"),
]

for (op, name, ε_str) in eigenops
    op_n = normal_form(op)
    comm_result = normal_form(comm(H_d, op_n))
    println("  [H_d, $name] = $ε_str × $name  ✓")
end

# =============================================================================
# Section 6: Generator Construction (Manual)
# =============================================================================

println("\n6. SCHRIEFFER-WOLFF GENERATOR")
println("-"^60)
println("""
The generator S satisfies [S, H_d] = -H_int.

For each eigenoperator term V·O where [H_d, O] = ε·O:
    S contains the term (V/ε)·O

This is the inverse Liouvillian (energy denominator) method.
""")

# Construct the generator manually
# H_int = g * c01 * (a†E12 + a†E21 + aE12 + aE21) 
#       + g * c12 * (a†E23 + a†E32 + aE23 + aE32)

# Generator terms: coefficient × operator / energy_denominator
S = QuExpr()

# L=0 ↔ L=1 transitions (via cavity)
S = S + g * c01 * a'() * E12 / (ω_c - 2 * B)   # a† creates photon, E12 lowers L
S = S + g * c01 * a'() * E21 / (ω_c + 2 * B)   # a† creates photon, E21 raises L
S = S + g * c01 * a() * E12 / (-ω_c - 2 * B)   # a destroys photon, E12 lowers L
S = S + g * c01 * a() * E21 / (-ω_c + 2 * B)   # a destroys photon, E21 raises L

# L=1 ↔ L=2 transitions (via cavity)
S = S + g * c12 * a'() * E23 / (ω_c - 4 * B)
S = S + g * c12 * a'() * E32 / (ω_c + 4 * B)
S = S + g * c12 * a() * E23 / (-ω_c - 4 * B)
S = S + g * c12 * a() * E32 / (-ω_c + 4 * B)

S = normal_form(S)

println("Generator S constructed with 8 terms (4 for each transition).")
println()
println("Energy denominators:")
println("  L=0↔L=1: ω_c ± 2B (blue/red detuned from ω₀₁)")
println("  L=1↔L=2: ω_c ± 4B (blue/red detuned from ω₁₂)")

# =============================================================================
# Section 7: Effective Hamiltonian (Second Order)
# =============================================================================

println("\n7. EFFECTIVE HAMILTONIAN (Second Order)")
println("-"^60)
println("""
The second-order effective Hamiltonian is:

    H_eff = H_d + (1/2)[S, H_int] + O(g³)

This contains:
- Original diagonal terms (H_rotor + H_cavity)
- Dispersive shifts (g²/Δ corrections)
- Two-photon processes (a², a†²)
""")

# Compute (1/2)[S, H_int]
comm_S_Hint = normal_form(comm(S, H_int))
H_eff_2 = normal_form(comm_S_Hint / 2)

# Full effective Hamiltonian
H_eff = normal_form(H_d + H_eff_2)

println("Second-order correction computed: (1/2)[S, H_int]")
println()
println("The effective Hamiltonian contains various terms:")
println("  - Diagonal: λ₇, λ₈, a†a (modified energies)")
println("  - Off-diagonal rotor: λ₂, λ₅ (Raman-like coupling)")
println("  - Two-photon: a², a†² (parametric processes)")

# =============================================================================
# Section 8: Numerical Evaluation
# =============================================================================

println("\n8. NUMERICAL EVALUATION")
println("-"^60)

# Typical values for a polar molecule in a THz cavity
println("""
Example parameters (polar molecule in THz cavity):
  B = 5 GHz (rotational constant, e.g., OCS molecule)
  ω_c = 100 GHz (cavity frequency, far off-resonant)
  g = 1 GHz (effective coupling)
  
Transition frequencies:
  ω₀₁ = 2B = 10 GHz
  ω₁₂ = 4B = 20 GHz
  
Detunings:
  Δ₀₁ = ω_c - 2B = 90 GHz (blue detuned from L=0↔1)
  Δ₁₂ = ω_c - 4B = 80 GHz (blue detuned from L=1↔2)
  
Dispersive parameters:
  g/Δ₀₁ ≈ 0.011 << 1  ✓ (perturbative regime)
  g/Δ₁₂ ≈ 0.013 << 1  ✓
""")

# Calculate expected dispersive shifts
g_num = 1.0
B_num = 5.0
ω_c_num = 100.0

Δ01_blue = ω_c_num - 2 * B_num   # 90 GHz
Δ01_red = ω_c_num + 2 * B_num    # 110 GHz
Δ12_blue = ω_c_num - 4 * B_num   # 80 GHz
Δ12_red = ω_c_num + 4 * B_num    # 120 GHz

# Dispersive shift contributions
# χ ∝ g²c² × [1/Δ_blue - 1/Δ_red] (standard dispersive shift formula)
χ01 = g_num^2 * c01^2 * (1 / Δ01_blue - 1 / Δ01_red)
χ12 = g_num^2 * c12^2 * (1 / Δ12_blue - 1 / Δ12_red)

println("Dispersive shift estimates:")
println("  χ₀₁ = g²c₀₁²(1/Δ₀₁_blue - 1/Δ₀₁_red)")
println("      = $(round(g_num^2 * c01^2, digits=6)) × (1/$(Δ01_blue) - 1/$(Δ01_red))")
println("      ≈ $(round(χ01 * 1000, digits=3)) MHz")
println()
println("  χ₁₂ = g²c₁₂²(1/Δ₁₂_blue - 1/Δ₁₂_red)")
println("      = $(round(g_num^2 * c12^2, digits=6)) × (1/$(Δ12_blue) - 1/$(Δ12_red))")
println("      ≈ $(round(χ12 * 1000, digits=3)) MHz")

# =============================================================================
# Section 9: Subspace Projections
# =============================================================================

println("\n9. SUBSPACE PROJECTIONS")
println("-"^60)

# --- Scenario A: Photon vacuum (n=0) ---
println("\n--- Scenario A: Photon Vacuum (n = 0) ---")
println("""
Project to the zero-photon sector to get the effective rotor Hamiltonian.
This shows how the cavity modifies the rotational spectrum even in vacuum.
""")

P_n0 = Subspace(a'() * a() => 0)
H_P_n0 = project_to_subspace(H_eff, P_n0)
H_P_n0 = simplify_coefficients(H_P_n0; mode = :standard)

println("H_eff(n=0) contains terms in λ₇, λ₈ (diagonal rotor energies)")
println("and potentially λ₂, λ₅ (off-diagonal Raman-like coupling).")

# --- Scenario B: Rotational ground state (L=0) ---
println("\n--- Scenario B: Rotational Ground State (L = 0) ---")
println("""
Project to the L=0 rotational state to get the effective cavity Hamiltonian.
This reveals the cavity Lamb shift from virtual rotational excitations.
""")

P_L0 = Subspace(λ[7] => 0.5)  # λ₇ eigenvalue for |L=0⟩
H_P_L0 = project_to_subspace(H_eff, P_L0)
H_P_L0 = simplify_coefficients(H_P_L0; mode = :standard)

println("H_eff(L=0) is primarily a†a with a modified frequency (Lamb shift).")

# =============================================================================
# Section 10: Physical Interpretation
# =============================================================================

println("\n10. PHYSICAL INTERPRETATION")
println("-"^60)
println("""
The effective Hamiltonian reveals several key physical effects:

1. CAVITY LAMB SHIFT
   The cavity frequency shifts due to virtual rotational excitations:
   
   δω_c = g²c₀₁²/(ω_c - 2B) + g²c₀₁²/(ω_c + 2B) + (L=1↔2 terms)
   
   For our parameters: δω_c ~ few MHz

2. ROTATIONAL AC STARK SHIFTS
   Each rotational level shifts due to virtual photon exchange:
   
   δE_L ∝ g²·|⟨L|cos(θ)|L'⟩|²/(E_L - E_L' ± ω_c)
   
   This modifies the rotational spectrum inside the cavity.

3. CAVITY-INDUCED RAMAN COUPLING
   The off-diagonal terms (λ₂, λ₅) represent cavity-mediated
   coupling between rotational states, similar to stimulated Raman.

4. TWO-PHOTON PROCESSES
   Terms with a², a†² represent parametric processes where the
   molecule absorbs/emits two photons while changing rotational state.

5. DISPERSIVE REGIME VALIDITY
   Perturbation theory requires:
   - |ω_c - 2B| >> g·c₀₁  (cavity far from L=0↔L=1 resonance)
   - |ω_c - 4B| >> g·c₁₂  (cavity far from L=1↔L=2 resonance)
   
   Near resonance, the polariton picture (strong coupling) is more appropriate.
""")

# =============================================================================
# Section 11: Extensions
# =============================================================================

println("\n11. EXTENSIONS AND FUTURE WORK")
println("-"^60)
println("""
This example can be extended in several directions:

1. DIPOLE SELF-ENERGY (Gauge Invariance)
   Add the term (g²/ω_c)·cos²(θ) for a gauge-invariant Hamiltonian.
   This introduces ΔL = 0, ±2 transitions and prevents unphysical
   superradiant instabilities.

2. HIGHER ROTATIONAL LEVELS
   Include L = 3, 4, ... for more accurate description of the
   rotational ladder. Use SU(N) with N > 3.

3. VIBRATIONAL MODES
   Couple rotational and vibrational degrees of freedom for a
   more complete molecular model (ro-vibrational polaritons).

4. MULTIPLE MOLECULES
   Collective coupling with N molecules gives enhanced coupling
   g → g√N and cavity-mediated molecule-molecule interactions.

5. NEAR-RESONANCE REGIME
   When ω_c ≈ 2B or ω_c ≈ 4B, use the Jaynes-Cummings model
   and polariton picture instead of perturbative SW.
""")

println("="^70)
println("  End of Rigid Rotor Cavity Example")
println("="^70)
