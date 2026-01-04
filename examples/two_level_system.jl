#=
Two-Level System with Transverse Field
=======================================

A simple two-level system with longitudinal (σz) and transverse (σx) fields:

    H = Δ/2 σz + ε σx

where:
- Δ: energy splitting (longitudinal field)
- ε: transverse field strength (perturbation)

When |ε| << |Δ|, we can use Schrieffer-Wolff to eliminate the off-diagonal
σx term, yielding an effective diagonal Hamiltonian with renormalized splitting:

    H_eff = (Δ/2 + ε²/Δ) σz + O(ε³)

This is equivalent to second-order perturbation theory, where the ground state
energy is lowered by -ε²/Δ and the excited state is raised by +ε²/Δ.

Exact solution: E_± = ±√(Δ²/4 + ε²) ≈ ±Δ/2 ± ε²/Δ for small ε

Reference: Standard quantum mechanics textbook (Sakurai, Griffiths, etc.)
=#

using UnitaryTransformations
using QuantumAlgebra
using Symbolics

println("="^60)
println("Two-Level System with Transverse Field")
println("="^60)

# Use σ± basis
QuantumAlgebra.use_σpm(true)

# Clear cached variables
UnitaryTransformations.clear_param_cache!()

# Define symbolic parameters
Δ = Pr"Δ"  # longitudinal splitting
ε = Pr"ε"  # transverse field (perturbation)

# Hamiltonian: H = Δ/2 σz + ε σx
# In σpm basis: σx = σ⁺ + σ⁻
H = Δ/2 * σz() + ε * (σp() + σm())

println("\n1. HAMILTONIAN")
println("-"^40)
println("H = Δ/2 σz + ε σx = Δ/2 σz + ε(σ⁺ + σ⁻)")
println("\nH = ", H)

# Define subspace: ground state |↓⟩ (σz = -1)
P = Subspace(σz() => -1)

println("\n2. SUBSPACE DEFINITION")
println("-"^40)
println("P = ground state subspace: σz → -1 (spin down)")

# Decompose
H_d, H_od = decompose(H, P)

println("\n3. HAMILTONIAN DECOMPOSITION")
println("-"^40)
println("H_diagonal     = ", H_d)
println("H_off-diagonal = ", H_od)

# Perform SW transformation
println("\n4. SCHRIEFFER-WOLFF TRANSFORMATION (order 2)")
println("-"^40)

result = schrieffer_wolff(H, P; order = 2)

println("Generator S = ", result.S)
println("\nEffective Hamiltonian H_eff = ", result.H_eff)

# Simplify and analyze
println("\n5. SIMPLIFIED H_eff")
println("-"^40)

H_eff_simp = simplify_coefficients(result.H_eff)
println("H_eff (simplified) = ", H_eff_simp)

terms = collect_terms(result.H_eff)
println("\nTerms in H_eff:")
for (op, coeff) in terms
    println("  ", coeff, "  ", op)
end

# Project to ground state
println("\n6. PROJECTED TO GROUND STATE")
println("-"^40)
println("H_P = ", result.H_P)

terms_P = collect_terms(result.H_P)
println("\nTerms in H_P:")
for (op, coeff) in terms_P
    println("  ", coeff, "  ", op)
end

# The ground state energy shift
println("\n7. ENERGY ANALYSIS")
println("-"^40)

# In the P subspace (σz = -1), we get the ground state energy
# E_g = -Δ/2 - ε²/Δ (lowered by second-order correction)
println("Ground state energy from H_P:")
for (op, coeff) in terms_P
    if op == "𝟙"
        println("  E_g = ", coeff)
    end
end

println("\nExpected from 2nd-order perturbation theory:")
println("  E_g = -Δ/2 - ε²/Δ")
println("  (Energy lowered by ε²/Δ due to mixing with excited state)")

# Exact solution comparison
println("\n8. COMPARISON WITH EXACT SOLUTION")
println("-"^40)
println("""
Exact eigenvalues: E_± = ±√(Δ²/4 + ε²)

Taylor expansion for small ε:
  E_± ≈ ±Δ/2 · √(1 + 4ε²/Δ²)
      ≈ ±Δ/2 · (1 + 2ε²/Δ²)
      = ±Δ/2 ± ε²/Δ

Ground state (E_-):
  E_g ≈ -Δ/2 - ε²/Δ  ✓ Matches our SW result!

Excited state (E_+):
  E_e ≈ +Δ/2 + ε²/Δ
""")

# Numerical verification
println("9. NUMERICAL VERIFICATION")
println("-"^40)
println("Parameters: Δ = 1.0, ε = 0.1 (perturbative regime: ε/Δ = 0.1)")

H_P_num = substitute_values(result.H_P, Dict(:Δ => 1.0, :ε => 0.1))
println("\nH_P(numerical) = ", H_P_num)

# Compute expected values
Δ_val = 1.0
ε_val = 0.1

E_g_exact = -sqrt(Δ_val^2/4 + ε_val^2)
E_g_SW = -Δ_val/2 - ε_val^2/Δ_val

println("\nGround state energy:")
println("  Exact:         E_g = ", E_g_exact)
println("  SW (2nd order): E_g = ", E_g_SW)
println(
    "  Error: ",
    abs(E_g_exact - E_g_SW),
    " (",
    100*abs(E_g_exact - E_g_SW)/abs(E_g_exact),
    "%)",
)

# Try with a larger perturbation to see breakdown
println("\n10. BREAKDOWN OF PERTURBATION THEORY")
println("-"^40)
println("Testing with larger ε to see when SW breaks down:")

for ε_test in [0.1, 0.3, 0.5, 0.8, 1.0]
    E_exact = -sqrt(Δ_val^2/4 + ε_test^2)
    E_SW = -Δ_val/2 - ε_test^2/Δ_val
    error_pct = 100*abs(E_exact - E_SW)/abs(E_exact)
    println("  ε/Δ = $(ε_test): Error = $(round(error_pct, digits=2))%")
end

println("\nAs expected, SW works well when ε << Δ and breaks down when ε ~ Δ")
