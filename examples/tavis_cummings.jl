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

This example uses the NEW SymSum type from QuantumAlgebra.jl which properly
handles symbolic sums with correct commutator semantics:
- Same-site terms: Σᵢ [Aᵢ, Bᵢ]
- Cross-site terms: Σᵢ Σⱼ≠ᵢ [Aᵢ, Bⱼ]

In the dispersive regime (|Δ| >> g√N), the Schrieffer-Wolff transformation
yields the effective Hamiltonian:

    H_eff = ω̃_c a†a + Σᵢ[ω̃_q/2 σz(i)] + χ·a†a·Jz + χ·Σᵢ≠ⱼ(σ⁺ᵢσ⁻ⱼ)

where:
- χ = g²/Δ is the single-atom dispersive shift
- Jz = Σᵢ σz(i)/2 is the collective spin operator
- The last term is the EXCHANGE INTERACTION between atoms

Key physics:
- Collective enhancement: The effective coupling scales as g√N
- Exchange interaction: Virtual photon exchange couples spin degrees of freedom
- Cavity-mediated entanglement: Foundation for multi-qubit gates

Reference: 
- Tavis & Cummings, Phys. Rev. 170, 379 (1968)
- Blais et al., PRA 75, 032329 (2007) for circuit QED context
=#

using QuantumAlgebra
using Symbolics

# The SymSum types are exported from QuantumAlgebra
# Import additional utilities for this example
import QuantumAlgebra: sumindex

println("="^70)
println("  Tavis-Cummings Model: N Atoms with Symbolic Sums")
println("="^70)

# Use σ± basis for cleaner algebra
QuantumAlgebra.use_σpm(true)

# Define symbolic parameters
@variables ω_c Δ g  # cavity freq, detuning (in rotating frame), coupling

println("\n1. HAMILTONIAN CONSTRUCTION WITH SYMBOLIC SUMS")
println("-"^60)
println("Using SymSum to represent Σᵢ with proper bound variable semantics")
println()

# Create a sum index
i = sumindex(1)

# Build the Tavis-Cummings Hamiltonian using SymSum
# H_cav = ω_c a†a (no sum - single cavity mode)
H_cav = ω_c * a'() * a()

# H_atom = Σᵢ (Δ/2) σz(i)
H_atom = SymSum(Δ/2 * σz(i), i)

# H_int = Σᵢ g(a†σ⁻(i) + aσ⁺(i))
H_int = SymSum(g * (a'()*σm(i) + a()*σp(i)), i)

# Total Hamiltonian (SymExpr combining QuExpr and SymSum terms)
H = SymExpr(H_cav) + H_atom + H_int

println("Cavity term:      H_cav = ω_c a†a")
println("Atom term:        H_atom = ", H_atom)
println("Interaction term: H_int  = ", H_int)
println()
println("Full Hamiltonian:")
println("  H = ", H)

# =============================================================================
# Section 2: Manual Schrieffer-Wolff for Symbolic Sums
# =============================================================================

println("\n2. SCHRIEFFER-WOLFF TRANSFORMATION (MANUAL)")
println("-"^60)
println("For the dispersive regime, we eliminate the atom-cavity coupling.")
println()

# The off-diagonal part is V = H_int
# The diagonal part is H_d = H_cav + H_atom
V = H_int

println("Off-diagonal perturbation:")
println("  V = ", V)

# The first-order generator S₁ satisfies [S₁, H_d] = -V
# For the JC/TC interaction: [S₁, H_d] = -g Σᵢ(a†σ⁻(i) + aσ⁺(i))
# 
# The solution is:
#   S₁ = Σᵢ (g/Δ)(aσ⁺(i) - a†σ⁻(i))
# 
# where we use [σz, σ±] = ±2σ± and [a†a, a] = -a, [a†a, a†] = a†

# Construct the generator
S1 = SymSum((g/Δ) * (a()*σp(i) - a'()*σm(i)), i)

println("\nFirst-order generator:")
println("  S₁ = ", S1)

# Verify: [S₁, H_d] should equal -V
# This requires computing commutators with SymSum
println("\nVerifying generator equation [S₁, H_d] = -V...")

# Commutator with cavity part: [S₁, ω_c a†a]
comm_cav = ω_c * comm(S1, a'()*a())
println("  [S₁, H_cav] = ", normal_form(comm_cav))

# Commutator with atom part: [S₁, Σⱼ (Δ/2)σz(j)]
# This is the key test! Uses same-site + cross-site decomposition
comm_atom = comm(S1, H_atom)
println("  [S₁, H_atom] = ", comm_atom)

# =============================================================================
# Section 3: Second-Order Correction - The Exchange Interaction
# =============================================================================

println("\n3. SECOND-ORDER CORRECTION: EXCHANGE INTERACTION")
println("-"^60)
println("The key physics: (1/2)[S₁, V] produces the exchange interaction!")
println()

# Compute [S₁, V]
# This is the crucial calculation that gives rise to the exchange term
comm_S1_V = comm(S1, V)

println("Commutator [S₁, V]:")
println("  ", comm_S1_V)

# Let's expand this for N=2 atoms to see the explicit terms
println("\n4. EXPANSION FOR N=2 ATOMS")
println("-"^60)

expanded_comm = expand_symbolic(comm_S1_V, 1:2)
expanded_normal = normal_form(expanded_comm)

println("Expanding [S₁, V] for 2 atoms:")
println("  ", expanded_normal)

# Simplify the coefficients
function simplify_expr(expr)
    result_terms = Dict{typeof(first(expr.terms)[1]),Any}()
    for (term, coeff) in expr.terms
        if coeff isa Symbolics.Num
            result_terms[term] = Symbolics.simplify(coeff)
        else
            result_terms[term] = coeff
        end
    end
    return QuantumAlgebra.QuExpr(result_terms)
end

expanded_simplified = simplify_expr(expanded_normal)
println("\nSimplified:")
println("  ", expanded_simplified)

# =============================================================================
# Section 5: Physical Interpretation
# =============================================================================

println("\n5. PHYSICAL INTERPRETATION")
println("-"^60)
println("""
The second-order effective Hamiltonian H_eff = H_d + (1/2)[S₁, V] contains:

1. CAVITY LAMB SHIFT
   The cavity frequency shifts due to virtual atom excitations.
   
2. AC STARK SHIFT  
   Each atom's frequency shifts by χ = g²/Δ due to virtual photon exchange.
   
3. DISPERSIVE COUPLING
   The term χ·a†a·Jz couples the photon number to the total spin.

4. EXCHANGE INTERACTION (NEW with SymSum!)
   The cross-site terms give:
   
     H_exchange = χ Σᵢ≠ⱼ (σ⁺ᵢσ⁻ⱼ + σ⁺ⱼσ⁻ᵢ)
   
   This represents:
   - Cavity-mediated spin-spin coupling
   - Virtual photon exchange: atom i emits, atom j absorbs
   - Foundation for cavity-mediated quantum gates
   - XY-type interaction: (σₓᵢσₓⱼ + σᵧᵢσᵧⱼ)/2

The exchange interaction was MISSING in the old approach because
sumindex(1) treated all sums as having the same index, giving only
same-site contributions.

With SymSum, we correctly get:
   [Σᵢ Aᵢ, Σⱼ Bⱼ] = Σᵢ[Aᵢ, Bᵢ] + Σᵢ Σⱼ≠ᵢ [Aᵢ, Bⱼ]
                     |_________|   |______________|
                      same-site     cross-site (exchange!)
""")

# =============================================================================
# Section 6: Comparison - Old vs New
# =============================================================================

println("\n6. COMPARISON: OLD (sumindex) vs NEW (SymSum)")
println("-"^60)

println("OLD APPROACH (broken for multi-atom physics):")
println("  i = sumindex(1)")
println("  V = g * (a'()*σm(i) + a()*σp(i))  # Just uses #₁")
println("  # Commutators treat all i's as the SAME site")
println("  # → Missing exchange terms!")
println()

println("NEW APPROACH (correct physics):")
println("  i = sumindex(1)")
println("  V = SymSum(g * (a'()*σm(i) + a()*σp(i)), i)")
println("  # SymSum tracks bound variable i")
println("  # Commutators split into same-site + cross-site")
println("  # → Exchange interaction emerges naturally!")

# =============================================================================
# Section 7: Explicit Two-Atom Calculation
# =============================================================================

println("\n7. EXPLICIT TWO-ATOM EFFECTIVE HAMILTONIAN")
println("-"^60)

# Build explicit two-atom Hamiltonian
@variables ω_c2 Δ2 g2

H2_cav = ω_c2 * a'() * a()
H2_atom = Δ2/2 * σz(1) + Δ2/2 * σz(2)
H2_int = g2 * (a'()*σm(1) + a()*σp(1)) + g2 * (a'()*σm(2) + a()*σp(2))

# Generator for two atoms
S2 = (g2/Δ2) * (a()*σp(1) - a'()*σm(1)) + (g2/Δ2) * (a()*σp(2) - a'()*σm(2))

# Compute (1/2)[S, V]
half_comm = (1//2) * normal_form(comm(S2, H2_int))

println("For 2 atoms, (1/2)[S₁, V] gives:")
for (term, coeff) in half_comm.terms
    coeff_simp = coeff isa Symbolics.Num ? Symbolics.simplify(coeff) : coeff
    println("  ", coeff_simp, " × ", term)
end

# Extract terms by type
println("\nGrouped by physics:")
println()

cavity_terms = []
single_atom_terms = []
exchange_terms = []

for (term, coeff) in half_comm.terms
    term_str = string(term)
    coeff_simp = coeff isa Symbolics.Num ? Symbolics.simplify(coeff) : coeff
    
    if occursin("a†", term_str) && occursin("a()", term_str)
        push!(cavity_terms, (term, coeff_simp))
    elseif occursin("σ⁺", term_str) && occursin("σ⁻", term_str)
        # Check if same site or different sites
        if (occursin("(1)", term_str) && occursin("(2)", term_str))
            push!(exchange_terms, (term, coeff_simp))
        else
            push!(single_atom_terms, (term, coeff_simp))
        end
    else
        push!(single_atom_terms, (term, coeff_simp))
    end
end

println("Cavity/dispersive terms:")
for (term, coeff) in cavity_terms
    println("  ", coeff, " × ", term)
end

println("\nSingle-atom terms:")
for (term, coeff) in single_atom_terms
    println("  ", coeff, " × ", term)
end

println("\nEXCHANGE TERMS (the key new physics!):")
for (term, coeff) in exchange_terms
    println("  ", coeff, " × ", term)
end

println("""

The exchange terms σ⁺(1)σ⁻(2) + σ⁺(2)σ⁻(1) represent:
  • Flip-flop interaction between atoms
  • Virtual photon exchange via the cavity
  • XY spin coupling: (σₓ₁σₓ₂ + σᵧ₁σᵧ₂)/2
  
This interaction enables:
  • iSWAP and √iSWAP gates
  • Entanglement generation
  • Collective spin dynamics
""")

# =============================================================================
# Section 8: Symbolic N-Atom Result
# =============================================================================

println("\n8. SYMBOLIC N-ATOM RESULT")
println("-"^60)

println("""
For N atoms, the effective Hamiltonian in the dispersive regime is:

  H_eff = ω̃_c a†a + Σᵢ ω̃_q/2 σz(i) + χ a†a Jz + χ Σᵢ≠ⱼ σ⁺ᵢσ⁻ⱼ

where:
  • χ = g²/Δ is the dispersive shift
  • Jz = (1/2)Σᵢ σz(i) is the collective spin z-component
  • ω̃_c = ω_c - Nχ/2 (collective Lamb shift)
  • ω̃_q = ω_q + χ (AC Stark shift)

In SymSum notation, the exchange term is:

  H_exchange = χ · Σ#₁(Σ#₂≠#₁(σ⁺(#₁)σ⁻(#₂)))

This nested sum with the j≠i constraint correctly represents:
  Σᵢ Σⱼ≠ᵢ σ⁺ᵢσ⁻ⱼ = Σᵢ<ⱼ (σ⁺ᵢσ⁻ⱼ + σ⁺ⱼσ⁻ᵢ)
""")

println("="^70)
println("  End of Tavis-Cummings Example")
println("="^70)
