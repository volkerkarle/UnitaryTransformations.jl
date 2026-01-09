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

# =============================================================================
# Section 1: System Setup
# =============================================================================

println("\n1. SYSTEM SETUP")
println("-"^60)

# 3-level rotor using transition operators
# States: |1⟩ = |L=0⟩, |2⟩ = |L=1⟩, |3⟩ = |L=2⟩
L = nlevel_ops(3, :L)

# Symbolic parameters
# Use symbolic coefficients for cleaner output
@variables B ω_c g c₀₁ c₁₂

# Rotor energies: E_L = B·L(L+1)
# E₀ = 0, E₁ = 2B, E₂ = 6B
H_rotor = 0 * L[1, 1] + 2B * L[2, 2] + 6B * L[3, 3]

# Cavity Hamiltonian
H_cav = ω_c * a'() * a()

# Dipole matrix elements (from Clebsch-Gordan coefficients)
# c₀₁ = 1/√3 ≈ 0.577   (⟨L=0|cos(θ)|L=1⟩)
# c₁₂ = √(2/15) ≈ 0.365 (⟨L=1|cos(θ)|L=2⟩)
# Keep these symbolic for cleaner algebra

# Dipole operator: cos(θ) with selection rule ΔL = ±1
cos_theta = c₀₁ * (L[1, 2] + L[2, 1]) + c₁₂ * (L[2, 3] + L[3, 2])

# Interaction: g·cos(θ)·(a + a†)
H_int = g * cos_theta * (a'() + a())

# Full Hamiltonian
H = normal_form(H_rotor + H_cav + H_int)

println("Rotor Hamiltonian: H_rotor = 0·|L=0⟩⟨L=0| + 2B·|L=1⟩⟨L=1| + 6B·|L=2⟩⟨L=2|")
println("Cavity Hamiltonian: H_cav = ω_c a†a")
println("Dipole operator: cos(θ) = c₀₁(|0⟩⟨1| + |1⟩⟨0|) + c₁₂(|1⟩⟨2| + |2⟩⟨1|)")
println("  c₀₁ = 1/√3 ≈ 0.577  (symbolic)")
println("  c₁₂ = √(2/15) ≈ 0.365  (symbolic)")
println("\nFull Hamiltonian:")
println("  H = ", H)

# =============================================================================
# Section 2: Schrieffer-Wolff Transformation
# =============================================================================

println("\n2. SCHRIEFFER-WOLFF TRANSFORMATION")
println("-"^60)

# --- Option A: Project to photon vacuum (n = 0) ---
println("\n--- Option A: Photon Vacuum (n = 0) ---")
println("Effective rotor Hamiltonian when cavity is in vacuum state.")
println("(Tip: set JULIA_NUM_THREADS>1 for faster order 4)")

P_n0 = Subspace(a'() * a() => 0)

# Order 2 (full SW)
t_o2 = @elapsed result_n0_o2 = schrieffer_wolff(H, P_n0; order = 2, simplify_mode = :fractions)
println("  Order 2 finished in $(round(t_o2, digits=2)) s")

# Order 4 using diagonal_only + parallel to save RAM/CPU
t_o4 = @elapsed result_n0_o4 = schrieffer_wolff(
    H,
    P_n0;
    order = 4,
    simplify_mode = :fast,
    diagonal_only = true,
    parallel = true,
)
println("  Order 4 (diagonal_only, parallel) finished in $(round(t_o4, digits=2)) s")

println("\nOrder 2 Generator S:")
println("  ", result_n0_o2.S)
println("\nOrder 2 Effective Hamiltonian H_eff:")
println("  ", result_n0_o2.H_eff)
println("\nOrder 2 Projected to n=0 subspace (H_P):")
println("  ", result_n0_o2.H_P)

println("\nOrder 4 Generator S (diagonal + off-diagonal up to g^4):")
println("  ", result_n0_o4.S)
println("\nOrder 4 Effective Hamiltonian H_eff:")
println("  ", result_n0_o4.H_eff)
println("\nOrder 4 Projected to n=0 subspace (H_P):")
println("  ", result_n0_o4.H_P)


# --- Option B: Project to rotational ground state (L = 0) ---
println("\n--- Option B: Rotational Ground State (L = 0) ---")
println("Effective cavity Hamiltonian when rotor is in ground state.")

P_L0 = Subspace(L[1, 1] => 1)

result_L0_o2 = schrieffer_wolff(H, P_L0; order = 2, simplify_mode = :fractions)
result_L0_o4 = schrieffer_wolff(H, P_L0; order = 4, simplify_mode = :fast, diagonal_only=true, parallel=true)

println("\nOrder 2 Generator S:")
println("  ", result_L0_o2.S)
println("\nOrder 2 Effective Hamiltonian H_eff:")
println("  ", result_L0_o2.H_eff)
println("\nOrder 2 Projected to L=0 subspace (H_P):")
println("  ", result_L0_o2.H_P)

println("\nOrder 4 Generator S (diagonal + off-diagonal up to g^4):")
println("  ", result_L0_o4.S)
println("\nOrder 4 Effective Hamiltonian H_eff:")
println("  ", result_L0_o4.H_eff)
println("\nOrder 4 Projected to L=0 subspace (H_P):")
println("  ", result_L0_o4.H_P)

# --- Option C: Project to both n=0 and L=0 (full ground state) ---
println("\n--- Option C: Full Ground State (n = 0, L = 0) ---")
println("Ground state energy shift from virtual excitations.")

P_ground = Subspace(a'() * a() => 0, L[1, 1] => 1)

result_ground = schrieffer_wolff(H, P_ground; order = 2, simplify_mode = :fractions)
result_L0n0_o4 = schrieffer_wolff(H, P_ground; order = 4, simplify_mode = :fast, diagonal_only=true, parallel=true)

println("\nGenerator S:")
println("  ", result_ground.S)
println("\nEffective Hamiltonian H_eff:")
println("  ", result_ground.H_eff)
println("\nProjected to ground state (H_P):")
println("  ", result_ground.H_P)

println("\nOrder 4 Generator S (diagonal + off-diagonal up to g^4):")
println("  ", result_L0n0_o4.S)
println("\nOrder 4 Effective Hamiltonian H_eff:")
println("  ", result_L0n0_o4.H_eff)
println("\nOrder 4 Projected to L=0 subspace (H_P):")
println("  ", result_L0n0_o4.H_P)

#


# =============================================================================
# Section 3: Simplification with Detunings
# =============================================================================

# Import types needed for coefficient manipulation
using QuantumAlgebra: QuExpr, QuTerm
import Symbolics.SymbolicUtils

println("\n3. SIMPLIFICATION WITH DETUNINGS")
println("-"^60)
println("Rewriting expressions in terms of physical detunings:")
println("  Δ₊ = ω_c + 2B  (blue detuning from L=0↔1)")
println("  Δ₋ = ω_c - 2B  (red detuning from L=0↔1)")
println("  Δ₁₂₊ = ω_c + 4B  (blue detuning from L=1↔2)")
println("  Δ₁₂₋ = ω_c - 4B  (red detuning from L=1↔2)")

# Define detuning symbols
@variables Δ₊ Δ₋ Δ₁₂₊ Δ₁₂₋

"""
    clean_rationals(expr)

Convert Rational{Int64} with denominator 1 to plain integers for cleaner display.
E.g., (12//1) → 12
"""
function clean_rationals(expr)
    # Use SymbolicUtils from Symbolics to traverse and clean up rationals
    uexpr = Symbolics.unwrap(expr)
    clean = SymbolicUtils.Postwalk(x -> (x isa Rational && isone(denominator(x))) ? numerator(x) : x; threaded=false)
    return Symbolics.wrap(clean(uexpr))
end

"""
    simplify_with_detunings(expr::QuExpr, B, ω_c, Δ₊, Δ₋, Δ₁₂₊, Δ₁₂₋)

Rewrite a QuExpr by substituting energy denominators with detuning symbols.

Substitution rules:
- ±2B ± ω_c patterns → ±Δ± (L=0↔1 transitions)
- ±4B ± ω_c patterns → ±Δ₁₂± (L=1↔2 transitions)
- Quadratic combinations like ω_c² - 4B² → Δ₊Δ₋
"""
function simplify_with_detunings(expr::QuExpr, B, ω_c, Δ₊, Δ₋, Δ₁₂₊, Δ₁₂₋)
    # Build substitution dictionary for linear detuning patterns
    # Note: We need to handle both expanded and factored forms
    linear_subs = Dict(
        # L=0↔1 transitions (±2B ± ω_c)
        ω_c + 2B => Δ₊,
        2B + ω_c => Δ₊,
        ω_c - 2B => Δ₋,
        -2B + ω_c => Δ₋,
        -ω_c - 2B => -Δ₊,
        -2B - ω_c => -Δ₊,
        -ω_c + 2B => -Δ₋,
        2B - ω_c => -Δ₋,
        # L=1↔2 transitions (±4B ± ω_c)
        ω_c + 4B => Δ₁₂₊,
        4B + ω_c => Δ₁₂₊,
        ω_c - 4B => Δ₁₂₋,
        -4B + ω_c => Δ₁₂₋,
        -ω_c - 4B => -Δ₁₂₊,
        -4B - ω_c => -Δ₁₂₊,
        -ω_c + 4B => -Δ₁₂₋,
        4B - ω_c => -Δ₁₂₋,
    )
    
    # Quadratic patterns (products of detunings)
    # ω_c² - 4B² = (ω_c + 2B)(ω_c - 2B) = Δ₊Δ₋
    # ω_c² - 16B² = (ω_c + 4B)(ω_c - 4B) = Δ₁₂₊Δ₁₂₋
    quadratic_subs = Dict(
        ω_c^2 - 4B^2 => Δ₊ * Δ₋,
        -4B^2 + ω_c^2 => Δ₊ * Δ₋,
        4B^2 - ω_c^2 => -Δ₊ * Δ₋,
        -ω_c^2 + 4B^2 => -Δ₊ * Δ₋,
        ω_c^2 - 16B^2 => Δ₁₂₊ * Δ₁₂₋,
        -16B^2 + ω_c^2 => Δ₁₂₊ * Δ₁₂₋,
        16B^2 - ω_c^2 => -Δ₁₂₊ * Δ₁₂₋,
        -ω_c^2 + 16B^2 => -Δ₁₂₊ * Δ₁₂₋,
        # Perfect square patterns from simplification.md
        # -16B² - 16Bω_c - 4ω_c² = -4(2B + ω_c)² = -4Δ₊²
        -16B^2 - 16B*ω_c - 4ω_c^2 => -4Δ₊^2,
        -4ω_c^2 - 16B*ω_c - 16B^2 => -4Δ₊^2,
        # -16B² + 16Bω_c - 4ω_c² = -4(ω_c - 2B)² = -4Δ₋²
        -16B^2 + 16B*ω_c - 4ω_c^2 => -4Δ₋^2,
        -4ω_c^2 + 16B*ω_c - 16B^2 => -4Δ₋^2,
        # -48B² - 48Bω_c - 12ω_c² = -12Δ₊²
        -48B^2 - 48B*ω_c - 12ω_c^2 => -12Δ₊^2,
        # -48B² + 48Bω_c - 12ω_c² = -12Δ₋²
        -48B^2 + 48B*ω_c - 12ω_c^2 => -12Δ₋^2,
        # -48B² + 12ω_c² = 12(ω_c² - 4B²) = 12Δ₊Δ₋
        -48B^2 + 12ω_c^2 => 12Δ₊ * Δ₋,
        12ω_c^2 - 48B^2 => 12Δ₊ * Δ₋,
        # Cross terms appearing at order 2:
        # 8B² + 2Bω_c - ω_c² = -(ω_c² - 2Bω_c - 8B²) = -(ω_c - 4B)(ω_c + 2B) = -Δ₁₂₋ Δ₊
        8B^2 + 2B*ω_c - ω_c^2 => -Δ₁₂₋ * Δ₊,
        -ω_c^2 + 2B*ω_c + 8B^2 => -Δ₁₂₋ * Δ₊,
        ω_c^2 - 2B*ω_c - 8B^2 => Δ₁₂₋ * Δ₊,
        # 8B² - 2Bω_c - ω_c² = -(ω_c² + 2Bω_c - 8B²) = -(ω_c + 4B)(ω_c - 2B) = -Δ₁₂₊ Δ₋
        8B^2 - 2B*ω_c - ω_c^2 => -Δ₁₂₊ * Δ₋,
        -ω_c^2 - 2B*ω_c + 8B^2 => -Δ₁₂₊ * Δ₋,
        ω_c^2 + 2B*ω_c - 8B^2 => Δ₁₂₊ * Δ₋,
        # -4B² - ω_c² + 4Bω_c = -(ω_c - 2B)² = -Δ₋² (already covered)
        # -4B² - ω_c² - 4Bω_c = -(ω_c + 2B)² = -Δ₊² (already covered)
    )
    
    # Process each term in the QuExpr
    result_terms = Dict{QuTerm,Any}()
    
    for (term, coeff) in expr.terms
        if coeff isa Num
            # Apply substitutions to the symbolic coefficient
            new_coeff = coeff
            
            # First expand to canonical form
            new_coeff = expand(new_coeff)
            
            # Apply quadratic substitutions first (they're more specific)
            for (pattern, replacement) in quadratic_subs
                new_coeff = Symbolics.substitute(new_coeff, pattern => replacement)
            end
            
            # Apply linear detuning substitutions
            for (pattern, replacement) in linear_subs
                new_coeff = Symbolics.substitute(new_coeff, pattern => replacement)
            end
            
            # Simplify fractions to combine terms
            new_coeff = simplify_fractions(new_coeff)
            
            # Clean up rationals like 12//1 → 12
            new_coeff = clean_rationals(new_coeff)
            
            result_terms[term] = new_coeff
        else
            result_terms[term] = coeff
        end
    end
    
    return QuExpr(result_terms)
end

# Apply simplification to all results
println("\n--- Simplified Option A: Photon Vacuum (n = 0) ---")
H_P_n0_simplified = simplify_with_detunings(result_n0_o4.H_P, B, ω_c, Δ₊, Δ₋, Δ₁₂₊, Δ₁₂₋)
println("Order 4 H_P (simplified):")
println("  ", H_P_n0_simplified)

println("\n--- Simplified Option B: Rotational Ground State (L = 0) ---")
H_P_L0_simplified = simplify_with_detunings(result_L0_o4.H_P, B, ω_c, Δ₊, Δ₋, Δ₁₂₊, Δ₁₂₋)
println("Order 4 H_P (simplified):")
println("  ", H_P_L0_simplified)

println("\n--- Simplified Option C: Full Ground State (n = 0, L = 0) ---")
H_P_ground_simplified = simplify_with_detunings(result_L0n0_o4.H_P, B, ω_c, Δ₊, Δ₋, Δ₁₂₊, Δ₁₂₋)
println("Order 4 H_P (simplified):")
println("  ", H_P_ground_simplified)

# Also show order-2 results simplified for comparison
println("\n--- Order 2 Results (simplified) ---")
H_P_n0_o2_simplified = simplify_with_detunings(result_n0_o2.H_P, B, ω_c, Δ₊, Δ₋, Δ₁₂₊, Δ₁₂₋)
println("Option A (n=0): ", H_P_n0_o2_simplified)

H_P_L0_o2_simplified = simplify_with_detunings(result_L0_o2.H_P, B, ω_c, Δ₊, Δ₋, Δ₁₂₊, Δ₁₂₋)
println("Option B (L=0): ", H_P_L0_o2_simplified)

H_P_ground_o2_simplified = simplify_with_detunings(result_ground.H_P, B, ω_c, Δ₊, Δ₋, Δ₁₂₊, Δ₁₂₋)
println("Option C (n=0,L=0): ", H_P_ground_o2_simplified)

# =============================================================================
# Section 4: Hermitian Polynomial Form
# =============================================================================

println("\n4. HERMITIAN POLYNOMIAL FORM")
println("-"^60)
println("Extracting coefficients for the effective cavity Hamiltonian:")
println("  H_eff = E₀ + A(a†² + a²) + Ω a†a + κ(a†⁴ + a⁴) + μ(a†³a + a†a³) + ν a†²a²")
println()

# Import extract_coefficient from the package
using UnitaryTransformations: extract_coefficient

"""
    extract_hermitian_form(H_eff::QuExpr)

Extract coefficients for the Hermitian polynomial form of an effective 
cavity Hamiltonian.

Returns a NamedTuple with:
- E0: constant (energy shift)
- A: coefficient of (a†² + a²) - squeezing parameter
- Ω: coefficient of a†a - effective frequency
- κ: coefficient of (a†⁴ + a⁴) - quartic squeezing
- μ: coefficient of (a†³a + a†a³) - three-photon process
- ν: coefficient of a†²a² - Kerr nonlinearity

Verifies Hermiticity by checking that conjugate pairs have equal coefficients.
"""
function extract_hermitian_form(H_eff::QuExpr)
    # Define the operator basis - use normal_form to ensure canonical ordering
    # Identity (constant term) is represented by empty operator product
    
    # Extract coefficients for each operator structure
    # For the constant term, we need to find terms with no operators
    const_coeff = nothing
    for (term, coeff) in H_eff.terms
        if isempty(term.bares.v)  # No bare operators = constant term
            const_coeff = coeff
            break
        end
    end
    
    # Quadratic operators
    coeff_a2 = extract_coefficient(H_eff, a()^2)           # a²
    coeff_adag2 = extract_coefficient(H_eff, a'()^2)       # a†²
    coeff_adaga = extract_coefficient(H_eff, a'()*a())     # a†a (number operator)
    
    # Quartic operators
    coeff_a4 = extract_coefficient(H_eff, a()^4)                    # a⁴
    coeff_adag4 = extract_coefficient(H_eff, a'()^4)                # a†⁴
    coeff_adag3a = extract_coefficient(H_eff, a'()^3*a())           # a†³a
    coeff_adaga3 = extract_coefficient(H_eff, a'()*a()^3)           # a†a³
    coeff_adag2a2 = extract_coefficient(H_eff, a'()^2*a()^2)        # a†²a²
    
    # Check Hermiticity: coefficients of conjugate pairs should be equal
    # (or both nothing)
    hermitian_ok = true
    warnings = String[]
    
    function check_equal(c1, c2, name1, name2)
        if c1 === nothing && c2 === nothing
            return true
        elseif c1 === nothing || c2 === nothing
            push!(warnings, "  ⚠ $name1 = $c1, $name2 = $c2 (one is missing)")
            return false
        else
            # For symbolic expressions, check structural equality
            if c1 isa Num && c2 isa Num
                diff = simplify(c1 - c2)
                if !isequal(diff, 0) && !isequal(Symbolics.unwrap(diff), 0)
                    push!(warnings, "  ⚠ $name1 - $name2 = $diff (should be 0)")
                    return false
                end
            elseif c1 != c2
                push!(warnings, "  ⚠ $name1 = $c1, $name2 = $c2 (not equal)")
                return false
            end
        end
        return true
    end
    
    hermitian_ok &= check_equal(coeff_adag2, coeff_a2, "coeff(a†²)", "coeff(a²)")
    hermitian_ok &= check_equal(coeff_adag4, coeff_a4, "coeff(a†⁴)", "coeff(a⁴)")
    hermitian_ok &= check_equal(coeff_adag3a, coeff_adaga3, "coeff(a†³a)", "coeff(a†a³)")
    
    if !hermitian_ok
        println("Hermiticity check warnings:")
        for w in warnings
            println(w)
        end
    end
    
    # Build the result - use the a†ⁿ coefficient as the canonical one for pairs
    return (
        E0 = const_coeff,
        A = coeff_adag2,      # squeezing: A(a†² + a²)
        Ω = coeff_adaga,      # frequency: Ω a†a  
        κ = coeff_adag4,      # quartic squeezing: κ(a†⁴ + a⁴)
        μ = coeff_adag3a,     # three-photon: μ(a†³a + a†a³)
        ν = coeff_adag2a2,    # Kerr: ν a†²a²
        hermitian = hermitian_ok,
    )
end

"""
    print_hermitian_form(coeffs::NamedTuple; name::String="H_eff")

Pretty-print the Hermitian polynomial form coefficients.
"""
function print_hermitian_form(coeffs::NamedTuple; name::String="H_eff")
    println("$name = E₀ + A(a†² + a²) + Ω a†a + κ(a†⁴ + a⁴) + μ(a†³a + a†a³) + ν a†²a²")
    println()
    
    function format_coeff(c)
        if c === nothing
            return "0"
        elseif c isa Num
            return string(clean_rationals(c))
        else
            return string(c)
        end
    end
    
    println("  E₀ = ", format_coeff(coeffs.E0))
    println("  A  = ", format_coeff(coeffs.A))
    println("  Ω  = ", format_coeff(coeffs.Ω))
    println("  κ  = ", format_coeff(coeffs.κ))
    println("  μ  = ", format_coeff(coeffs.μ))
    println("  ν  = ", format_coeff(coeffs.ν))
    println()
    println("  Hermitian: ", coeffs.hermitian ? "✓" : "✗")
end

"""
    extract_rotor_hermitian_form(H_eff::QuExpr, L)

Extract coefficients for the Hermitian polynomial form of an effective 
rotor Hamiltonian (3-level system: L=0,1,2).

Returns a NamedTuple with diagonal and off-diagonal terms.
"""
function extract_rotor_hermitian_form(H_eff::QuExpr, L)
    # Diagonal terms (populations/energies)
    E0 = extract_coefficient(H_eff, L[1,1])   # |L=0⟩⟨L=0|
    E1 = extract_coefficient(H_eff, L[2,2])   # |L=1⟩⟨L=1|
    E2 = extract_coefficient(H_eff, L[3,3])   # |L=2⟩⟨L=2|
    
    # Off-diagonal terms (coherences/couplings)
    # L=0 ↔ L=1 coupling
    V01 = extract_coefficient(H_eff, L[1,2])  # |L=0⟩⟨L=1|
    V10 = extract_coefficient(H_eff, L[2,1])  # |L=1⟩⟨L=0|
    
    # L=1 ↔ L=2 coupling  
    V12 = extract_coefficient(H_eff, L[2,3])  # |L=1⟩⟨L=2|
    V21 = extract_coefficient(H_eff, L[3,2])  # |L=2⟩⟨L=1|
    
    # L=0 ↔ L=2 coupling (second-order process, ΔL=2)
    V02 = extract_coefficient(H_eff, L[1,3])  # |L=0⟩⟨L=2|
    V20 = extract_coefficient(H_eff, L[3,1])  # |L=2⟩⟨L=0|
    
    return (
        E0 = E0, E1 = E1, E2 = E2,
        V01 = V01, V10 = V10,
        V12 = V12, V21 = V21,
        V02 = V02, V20 = V20,
    )
end

"""
    print_rotor_hermitian_form(coeffs::NamedTuple; name::String="H_eff")

Pretty-print the rotor Hermitian polynomial form coefficients.
"""
function print_rotor_hermitian_form(coeffs::NamedTuple; name::String="H_eff")
    println("$name = E₀|0⟩⟨0| + E₁|1⟩⟨1| + E₂|2⟩⟨2| + V₀₁(|0⟩⟨1| + h.c.) + V₁₂(|1⟩⟨2| + h.c.) + V₀₂(|0⟩⟨2| + h.c.)")
    println()
    
    function format_coeff(c)
        if c === nothing
            return "0"
        elseif c isa Num
            return string(clean_rationals(c))
        else
            return string(c)
        end
    end
    
    println("Diagonal (energies):")
    println("  E₀ = ", format_coeff(coeffs.E0))
    println("  E₁ = ", format_coeff(coeffs.E1))
    println("  E₂ = ", format_coeff(coeffs.E2))
    println()
    println("Off-diagonal (couplings):")
    println("  V₀₁ = ", format_coeff(coeffs.V01), "  (ΔL=1)")
    println("  V₁₂ = ", format_coeff(coeffs.V12), "  (ΔL=1)")
    println("  V₀₂ = ", format_coeff(coeffs.V02), "  (ΔL=2, cavity-mediated)")
end

# =============================================================================
# Apply to Option A (n=0, photon vacuum) - Effective ROTOR Hamiltonian
# =============================================================================

println("--- Option A: Photon Vacuum (n = 0) ---")
println("Effective rotor Hamiltonian when cavity is in vacuum state.")
println()

println("Order 4 H_P (simplified with detunings):")
coeffs_n0_simplified = extract_rotor_hermitian_form(H_P_n0_simplified, L)
print_rotor_hermitian_form(coeffs_n0_simplified; name="H_eff")

println()
println("Order 2 H_P (simplified with detunings):")
coeffs_n0_o2 = extract_rotor_hermitian_form(H_P_n0_o2_simplified, L)
print_rotor_hermitian_form(coeffs_n0_o2; name="H_eff (order 2)")

# =============================================================================
# Apply to Option B (L=0, rotational ground state) - Effective CAVITY Hamiltonian
# =============================================================================

println()
println("--- Option B: Rotational Ground State (L = 0) ---")
println("Effective cavity Hamiltonian when rotor is in ground state.")
println()

# First extract from the simplified version
println("From simplified H_P (with detunings):")
coeffs_L0_simplified = extract_hermitian_form(H_P_L0_simplified)
print_hermitian_form(coeffs_L0_simplified; name="H_eff")

# Also show the order-2 result for comparison
println()
println("--- Order 2 Comparison ---")
coeffs_L0_o2 = extract_hermitian_form(H_P_L0_o2_simplified)
print_hermitian_form(coeffs_L0_o2; name="H_eff (order 2)")

# =============================================================================
# Section 5: Physical Interpretation
# =============================================================================

println("\n5. PHYSICAL INTERPRETATION")
println("-"^60)
println("""
The effective cavity Hamiltonian for a rigid rotor in its ground state (L=0):

  H_eff = E₀ + A(a†² + a²) + Ω a†a + κ(a†⁴ + a⁴) + μ(a†³a + a†a³) + ν a†²a²

Physical meaning of each term:

• E₀: Lamb shift - ground state energy shift from virtual excitations

• A(a†² + a²): SQUEEZING - the cavity ground state becomes a squeezed state!
  This arises because the rotor-cavity coupling breaks time-reversal symmetry
  at order g². Non-zero when Δ₊ ≠ Δ₋ (i.e., ω_c ≠ 0).

• Ω a†a: Frequency renormalization - AC Stark shift of the cavity.
  Contains both dispersive shift and Bloch-Siegert-like corrections.

• κ(a†⁴ + a⁴): Quartic squeezing - higher-order ground state deformation.

• μ(a†³a + a†a³): Three-photon nonlinearity - couples different Fock states.

• ν a†²a²: Kerr nonlinearity - photon-photon interaction.
  Leads to photon blockade and non-classical light generation.

Key insight: The factors (Δ₊ - Δ₋) = 4B appear throughout, reflecting the
asymmetry between red and blue sidebands due to the rotational structure.
""")

println("="^70)
println("  End of Rigid Rotor Cavity Example")
println("="^70)
