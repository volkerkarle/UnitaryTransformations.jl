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

Hamiltonian:
    H = E₁|1⟩⟨1| + E₂|2⟩⟨2| + E₃|3⟩⟨3| + g₁(|1⟩⟨3| + |3⟩⟨1|) + g₂(|2⟩⟨3| + |3⟩⟨2|)

where:
- E₁, E₂: Ground state energies (can be degenerate or split)
- E₃: Excited state energy (large detuning Δ = E₃ - E₁,₂)
- g₁, g₂: Coupling strengths to the excited state

In the limit Δ >> g₁, g₂, the Schrieffer-Wolff transformation eliminates
the excited state, giving an effective two-level system with:
- AC Stark shifts of ground states
- Two-photon Raman coupling ∝ g₁g₂/Δ between |1⟩ and |2⟩

Physical applications:
- Raman transitions in laser cooling
- Electromagnetically induced transparency (EIT)
- Two-photon processes in NMR
- Lambda systems in quantum optics

Reference: Cohen-Tannoudji, Arimondo, "Laser cooling and trapping of neutral atoms"
=#

using UnitaryTransformations
using QuantumAlgebra
using Symbolics

println("="^60)
println("Three-Level Atom (Lambda Configuration)")
println("="^60)

# =============================================================================
# Section 1: System Setup
# =============================================================================

println("\n1. SYSTEM SETUP")
println("-"^40)

# 3-level atom using transition operators
# States: |1⟩, |2⟩ = ground states, |3⟩ = excited state
σ = nlevel_ops(3, :σ)

# Symbolic parameters
@variables E₁ E₂ E₃ g₁ g₂

# Atomic Hamiltonian (diagonal)
H_atom = E₁ * σ[1, 1] + E₂ * σ[2, 2] + E₃ * σ[3, 3]

# Coupling to excited state (off-diagonal)
H_coupling = g₁ * (σ[1, 3] + σ[3, 1]) + g₂ * (σ[2, 3] + σ[3, 2])

# Full Hamiltonian
H = normal_form(H_atom + H_coupling)

println("Three-level Lambda system:")
println("  |3⟩ = excited state (energy E₃)")
println("  |1⟩, |2⟩ = ground states (energies E₁, E₂)")
println()
println("H_atom = E₁|1⟩⟨1| + E₂|2⟩⟨2| + E₃|3⟩⟨3|")
println("H_coupling = g₁(|1⟩⟨3| + h.c.) + g₂(|2⟩⟨3| + h.c.)")
println()
println("Full Hamiltonian:")
println("  H = ", H)

# =============================================================================
# Section 2: Schrieffer-Wolff Transformation
# =============================================================================

println("\n2. SCHRIEFFER-WOLFF TRANSFORMATION")
println("-"^40)

# --- Option A: Project out excited state |3⟩ ---
println("\n--- Option A: Eliminate Excited State ---")
println("Project to ground state manifold {|1⟩, |2⟩} by eliminating |3⟩.")

# The excited state has σ[3,3] = 0 in the ground manifold
P_ground = Subspace(σ[3, 3] => 0)

result_ground = schrieffer_wolff(H, P_ground; order = 2, simplify_mode = :fractions)

println("\nGenerator S:")
println("  ", result_ground.S)
println("\nEffective Hamiltonian H_eff:")
println("  ", result_ground.H_eff)
println("\nProjected to ground manifold (H_P):")
println("  ", result_ground.H_P)

# --- Option B: Project to state |1⟩ only ---
println("\n--- Option B: Project to State |1⟩ ---")
println("Get the energy of state |1⟩ including all virtual corrections.")

P_1 = Subspace(σ[1, 1] => 1)

result_1 = schrieffer_wolff(H, P_1; order = 2, simplify_mode = :fractions)

println("\nProjected to |1⟩ (H_P):")
println("  ", result_1.H_P)

# --- Option C: Project to state |2⟩ only ---
println("\n--- Option C: Project to State |2⟩ ---")
println("Get the energy of state |2⟩ including all virtual corrections.")

P_2 = Subspace(σ[2, 2] => 1)

result_2 = schrieffer_wolff(H, P_2; order = 2, simplify_mode = :fractions)

println("\nProjected to |2⟩ (H_P):")
println("  ", result_2.H_P)

# =============================================================================
# Section 3: Numerical Evaluation
# =============================================================================

println("\n3. NUMERICAL EVALUATION")
println("-"^40)

println("""
Example parameters (typical atomic system):
  E₁ = 0 (reference)
  E₂ = 0.1 (small ground state splitting)
  E₃ = 10 (large detuning from excited state)
  g₁ = 0.5, g₂ = 0.3 (coupling strengths)
  
Detunings:
  Δ₁ = E₃ - E₁ = 10
  Δ₂ = E₃ - E₂ = 9.9
  
Perturbative parameter:
  g/Δ ~ 0.05 << 1  (valid regime)
""")

# Substitute numerical values
params = Dict(:E₁ => 0.0, :E₂ => 0.1, :E₃ => 10.0, :g₁ => 0.5, :g₂ => 0.3)

H_P_ground_num = substitute_values(result_ground.H_P, params)
H_P_1_num = substitute_values(result_1.H_P, params)
H_P_2_num = substitute_values(result_2.H_P, params)

println("Numerical results:")
println("\n  Option A (ground manifold):")
println("    H_P = ", H_P_ground_num)
println("\n  Option B (|1⟩): H_P = ", H_P_1_num)
println("\n  Option C (|2⟩): H_P = ", H_P_2_num)

# =============================================================================
# Section 4: Physical Interpretation
# =============================================================================

println("\n4. PHYSICAL INTERPRETATION")
println("-"^40)
println("""
The effective Hamiltonian in the ground manifold contains:

1. AC STARK SHIFTS
   Each ground state |i⟩ shifts by -gᵢ²/(E₃ - Eᵢ)
   This is the light shift from virtual excitation to |3⟩.

2. RAMAN COUPLING
   An effective coupling appears between |1⟩ and |2⟩:
   
     H_Raman ∝ g₁g₂/Δ × (|1⟩⟨2| + |2⟩⟨1|)
   
   This is a TWO-PHOTON process: |1⟩ → |3⟩ → |2⟩

3. EFFECTIVE TWO-LEVEL SYSTEM
   After eliminating |3⟩, we have an effective qubit in {|1⟩, |2⟩}
   with modified energies and Raman coupling.

Applications:
- Raman cooling: Drive |1⟩ ↔ |2⟩ transitions via Raman
- EIT: Interference between excitation pathways
- Quantum gates: Use Raman coupling for qubit rotations
""")

# =============================================================================
# Section 5: Verify Raman Coupling Term
# =============================================================================

println("\n5. RAMAN COUPLING ANALYSIS")
println("-"^40)

# Extract the |1⟩⟨2| coefficient from the effective Hamiltonian
H_P = result_ground.H_P
raman_coeff = extract_coefficient(H_P, σ[1, 2])

println("Coefficient of |1⟩⟨2| (Raman coupling):")
println("  ", raman_coeff)
println()
println("Expected form: g₁g₂ × f(E₁, E₂, E₃)")
println("where f contains energy denominators from virtual |3⟩ excitation.")

println("\n" * "="^60)
println("End of Three-Level Atom Example")
println("="^60)
