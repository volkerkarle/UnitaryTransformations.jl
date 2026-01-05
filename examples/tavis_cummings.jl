#=
Tavis-Cummings Model: N Atoms in a Cavity
==========================================

The Tavis-Cummings Hamiltonian generalizes the Jaynes-Cummings model to N 
identical two-level atoms coupled to a single cavity mode:

    H = ω_c a†a + Σᵢ[ω_q/2 σz(i)] + g·Σᵢ[a†σ⁻(i) + a σ⁺(i)]

where:
- ω_c: cavity frequency
- ω_q: qubit frequency (identical for all atoms)
- g: single-atom coupling strength
- Δ = ω_q - ω_c: detuning
- i ∈ {1, 2, ..., N}: atom index

In QuantumAlgebra, the sum over atoms is represented using a "sum index" #₁,
which symbolically represents Σᵢ.

In the dispersive regime (|Δ| >> g√N), the Schrieffer-Wolff transformation
yields the effective Hamiltonian:

    H_eff = ω̃_c a†a + Σᵢ[ω̃_q/2 σz(i)] + χ·a†a·Jz

where:
- χ = g²/Δ is the single-atom dispersive shift
- Jz = Σᵢ σz(i)/2 is the collective spin operator
- ω̃_c = ω_c - Nχ (collective Lamb shift)
- ω̃_q = ω_q + 2χ (AC Stark shift)

Key physics:
- Collective enhancement: The effective coupling scales as g√N
- Superradiance: Collective states decay at enhanced rates
- Collective dispersive shift: Cavity shift depends on total spin

Reference: 
- Tavis & Cummings, Phys. Rev. 170, 379 (1968)
- Blais et al., PRA 75, 032329 (2007) for circuit QED context
=#

using UnitaryTransformations
using QuantumAlgebra
using Symbolics

println("="^65)
println("  Tavis-Cummings Model: N Atoms in a Cavity")
println("="^65)

# Use σ± basis for cleaner algebra
QuantumAlgebra.use_σpm(true)

# Clear any cached symbolic variables from previous runs
UnitaryTransformations.clear_param_cache!()

# Create a sum index to represent Σᵢ over N atoms
# In QuantumAlgebra, sumindex(1) creates the symbolic index #₁
i = QuantumAlgebra.sumindex(1)

# Define symbolic parameters
@variables ω_c Δ g  # cavity freq, detuning, coupling

println("\n1. HAMILTONIAN")
println("-"^50)
println("H = ω_c a†a + Σᵢ[Δ/2 σz(i)] + g·Σᵢ[a†σ⁻(i) + a σ⁺(i)]")
println()

# Build the Tavis-Cummings Hamiltonian
# The sum index #₁ automatically represents the sum over all atoms
H_cav = ω_c * a'() * a()           # Cavity energy
H_atom = Δ/2 * σz(i)               # Atom energies (summed over i)
H_int = g * (a'()*σm(i) + a()*σp(i))  # Jaynes-Cummings coupling (summed)

H = H_cav + H_atom + H_int

println("In QuantumAlgebra notation:")
println("  H = ", normal_form(H))
println()
println("Note: #₁ represents the sum index Σᵢ over N atoms")

# Define subspace: cavity in vacuum state (dispersive regime)
P = Subspace(a'()*a() => 0)

println("\n2. SUBSPACE DEFINITION")
println("-"^50)
println("P = {|n⟩ : n = 0} (cavity vacuum / low photon number)")
println()
println("This is appropriate for the dispersive regime where")
println("photon number is small and we want to eliminate the")
println("direct atom-cavity coupling.")

# Decompose into diagonal and off-diagonal parts
H_d, H_od = decompose(H, P)

println("\n3. HAMILTONIAN DECOMPOSITION")
println("-"^50)
println("H_diagonal     = ", H_d)
println("H_off-diagonal = ", H_od)

# Perform Schrieffer-Wolff transformation to second order
println("\n4. SCHRIEFFER-WOLFF TRANSFORMATION")
println("-"^50)

result = schrieffer_wolff(H, P; order = 2)

println("Generator:")
S_simplified = simplify_coefficients(result.S)
println("  S = ", S_simplified)

println("\nEffective Hamiltonian (order 2):")
H_eff = simplify_coefficients(result.H_eff)
println("  H_eff = ", H_eff)

# Verify the first-order generator equation
println("\n5. VERIFICATION")
println("-"^50)
# Get first-order generator
S1 = solve_for_generator(H_d, H_od, P)
residual = normal_form(comm(S1, H_d) + H_od)
# Each coefficient should simplify to zero
verified = isempty(residual.terms)
if !verified
    # Check if coefficients simplify to zero
    verified = all(iszero(Symbolics.simplify(c)) for (_, c) in residual.terms)
end
if verified
    println("First-order generator equation [S₁, H_d] = -H_od: ✓ Verified")
    println("  S₁ = ", S1)
else
    println("Generator equation not satisfied. Residual: ", residual)
end

# Physical interpretation
println("\n6. PHYSICAL INTERPRETATION")
println("-"^50)
println("""
The effective Hamiltonian terms can be written in standard form:

  H_eff = (ω_c - χ·N)a†a + Σᵢ[(Δ + 2χ)/2 σz(i)] + 2χ·a†a·Σᵢ[σ⁺(i)σ⁻(i)]

where χ = g²/(Δ - ω_c) ≈ g²/Δ is the dispersive shift.

Converting to Jz = Σᵢ σz(i)/2 (collective spin):

  H_eff = ω̃_c·a†a + ω̃_q·Jz + χ·a†a·Jz

Key physics of the Tavis-Cummings model:

1. COLLECTIVE LAMB SHIFT
   The cavity frequency shifts by -χ per atom (total -Nχ).
   This is a collective effect where all atoms contribute.

2. AC STARK SHIFT  
   Each atom's frequency shifts by 2χ due to virtual photon exchange.

3. DISPERSIVE COUPLING
   The term χ·a†a·Jz couples the photon number to the total spin.
   - Cavity frequency depends on collective spin state
   - Individual atom frequencies shift with photon number

4. COLLECTIVE ENHANCEMENT
   For N atoms, the vacuum Rabi splitting scales as g√N,
   leading to enhanced dispersive effects.

This is the foundation for:
- Multi-qubit dispersive readout in circuit QED
- Collective spin squeezing
- Quantum non-demolition measurements
""")

# Numerical example
println("\n7. NUMERICAL EXAMPLE")
println("-"^50)
println("Parameters: Δ = 1.0 GHz, g = 0.1 GHz (dispersive: g/Δ = 0.1)")
println()
println("Dispersive shift χ = g²/Δ = 0.01 GHz = 10 MHz")
println()
println("For N = 10 atoms:")
println("  - Collective Lamb shift: 10 × 10 MHz = 100 MHz")
println("  - AC Stark shift per atom: 20 MHz")
println("  - Dispersive coupling: 10 MHz × a†a × Jz")
