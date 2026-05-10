#=
Rigid Rotor in a Cavity (Detuning Basis, Option B)
==================================================

Runs the L=0 projection in a detuning parametrization from the start:

    Δ₊ = ω_c + 2B
    Δ₋ = ω_c - 2B

so the SW expansion is carried out directly in terms of Δ₊, Δ₋.

This keeps expressions smaller and avoids the large post-SW detuning rewrite.
=#

using UnitaryTransformations
using QuantumAlgebra
using Symbolics

import Symbolics.SymbolicUtils
using UnitaryTransformations: extract_coefficient

println("="^70)
println("  Rigid Rotor in a Cavity: SW in Detuning Basis")
println("="^70)
println("(Tip: set JULIA_NUM_THREADS to use all cores)")

PRINT_FULL_HP = get(ENV, "UT_PRINT_FULL_HP", "0") == "1"
PARALLEL_SW = get(ENV, "UT_PARALLEL_SW", "0") == "1"
G_SIGN = parse(Float64, get(ENV, "UT_G_SIGN", "1"))

# =============================================================================
# Detuning parametrization
# =============================================================================

@variables Δ₊ Δ₋ g c₀₁ c₁₂

# Express ω_c and B in terms of detunings
ω_c = (Δ₊ + Δ₋) / 2
B = (Δ₊ - Δ₋) / 4

# Secondary detunings for L=1↔2 transitions (derived)
Δ12p_expr = (3Δ₊ - Δ₋) / 2
Δ12m_expr = (-Δ₊ + 3Δ₋) / 2
@variables Δ₁₂₊ Δ₁₂₋

# =============================================================================
# System setup
# =============================================================================

L = nlevel_ops(3, :L)

H_rotor = 0 * L[1, 1] + 2B * L[2, 2] + 6B * L[3, 3]
H_cav = ω_c * a'() * a()

cos_theta = c₀₁ * (L[1, 2] + L[2, 1]) + c₁₂ * (L[2, 3] + L[3, 2])
H_int = (G_SIGN * g) * cos_theta * (a'() + a())

H = normal_form(H_rotor + H_cav + H_int)

println("Detunings: Δ₊ = ω_c + 2B, Δ₋ = ω_c - 2B")
println("Derived:   Δ₁₂₊ = ω_c + 4B = (3Δ₊ - Δ₋)/2")
println("           Δ₁₂₋ = ω_c - 4B = (-Δ₊ + 3Δ₋)/2")

# =============================================================================
# SW transformation (Option B: L=0)
# =============================================================================

P_L0 = Subspace(L[1, 1] => 1)

t_o4 = @elapsed result_L0_o4 = schrieffer_wolff(
    H,
    P_L0;
    order = 4,
    simplify_mode = :fast,
    parallel = PARALLEL_SW,
)

println("Order 4 finished in $(round(t_o4, digits=2)) s (parallel=$(PARALLEL_SW))")
println("\nOrder 4 Projected to L=0 subspace (H_P):")
println("  terms = ", length(result_L0_o4.H_P.terms))
if PRINT_FULL_HP
    println("  ", result_L0_o4.H_P)
end

# =============================================================================
# Hermitian polynomial form extraction
# =============================================================================

function clean_rationals(expr)
    uexpr = Symbolics.unwrap(expr)
    clean = SymbolicUtils.Postwalk(
        x -> (x isa Rational && isone(denominator(x))) ? numerator(x) : x;
        threaded = false,
    )
    return Symbolics.wrap(clean(uexpr))
end

function extract_hermitian_form(H_eff::QuExpr)
    const_coeff = nothing
    for (term, coeff) in H_eff.terms
        if isempty(term.bares.v)
            const_coeff = coeff
            break
        end
    end

    coeff_a2 = extract_coefficient(H_eff, a()^2)
    coeff_adag2 = extract_coefficient(H_eff, a'()^2)
    coeff_adaga = extract_coefficient(H_eff, a'() * a())

    coeff_a4 = extract_coefficient(H_eff, a()^4)
    coeff_adag4 = extract_coefficient(H_eff, a'()^4)
    coeff_adag3a = extract_coefficient(H_eff, a'()^3 * a())
    coeff_adaga3 = extract_coefficient(H_eff, a'() * a()^3)
    coeff_adag2a2 = extract_coefficient(H_eff, a'()^2 * a()^2)

    hermitian_ok = true
    warnings = String[]

    function numeric_zero(expr)
        subs_list = [
            Dict(Δ₊ => 7.0, Δ₋ => 5.0, g => 0.1, c₀₁ => 0.6, c₁₂ => 0.4),
            Dict(Δ₊ => 9.0, Δ₋ => 4.0, g => 0.2, c₀₁ => 0.5, c₁₂ => 0.3),
        ]
        for subs in subs_list
            val = Symbolics.substitute(expr, subs)
            val = val isa Num ? Symbolics.value(val) : val
            if val isa Num
                return false
            end
            if !isfinite(val) || abs(val) > 1e-8
                return false
            end
        end
        return true
    end

    function check_equal(c1, c2, name1, name2)
        if c1 === nothing && c2 === nothing
            return true
        elseif c1 === nothing || c2 === nothing
            push!(warnings, "  ⚠ $name1 = $c1, $name2 = $c2 (one is missing)")
            return false
        else
            if c1 isa Num && c2 isa Num
                diff = simplify_fractions(expand(c1 - c2))
                diff = clean_rationals(diff)
                if !isequal(diff, 0) && !isequal(Symbolics.unwrap(diff), 0)
                    if !numeric_zero(diff)
                        push!(warnings, "  ⚠ $name1 - $name2 = $diff (should be 0)")
                        return false
                    end
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

    return (
        E0 = const_coeff,
        A = coeff_adag2,
        Ω = coeff_adaga,
        κ = coeff_adag4,
        μ = coeff_adag3a,
        ν = coeff_adag2a2,
        hermitian = hermitian_ok,
    )
end

function rewrite_secondary_detunings(coeff)
    coeff === nothing && return nothing
    coeff isa Num || return coeff
    expanded = expand(coeff)
    rewritten = Symbolics.substitute(
        expanded,
        Dict(Δ12p_expr => Δ₁₂₊, Δ12m_expr => Δ₁₂₋),
    )
    return clean_rationals(rewritten)
end

function rewrite_coeffs_secondary_detunings(coeffs::NamedTuple)
    return (
        E0 = rewrite_secondary_detunings(coeffs.E0),
        A = rewrite_secondary_detunings(coeffs.A),
        Ω = rewrite_secondary_detunings(coeffs.Ω),
        κ = rewrite_secondary_detunings(coeffs.κ),
        μ = rewrite_secondary_detunings(coeffs.μ),
        ν = rewrite_secondary_detunings(coeffs.ν),
        hermitian = coeffs.hermitian,
    )
end

function print_hermitian_form(coeffs::NamedTuple; name::String = "H_eff")
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

println("\n--- Hermitian Polynomial Form (Option B, detuning basis) ---")
t_coeffs = @elapsed coeffs_L0 = extract_hermitian_form(result_L0_o4.H_P)
println("  Coefficient extraction in $(round(t_coeffs, digits=2)) s")
t_rewrite = @elapsed coeffs_L0 = rewrite_coeffs_secondary_detunings(coeffs_L0)
println("  Secondary detuning rewrite in $(round(t_rewrite, digits=2)) s")
print_hermitian_form(coeffs_L0; name = "H_eff")
