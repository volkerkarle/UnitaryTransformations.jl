#=
Compare include_QQ on/off for rigid rotor (Option B)
====================================================

Runs SW4 twice (include_QQ=false/true), extracts the Hermitian coefficients,
and checks whether the differences vanish when c12 = 0.
=#

using UnitaryTransformations
using QuantumAlgebra
using Symbolics

import Symbolics.SymbolicUtils
using UnitaryTransformations: extract_coefficient

println("="^70)
println("  Compare include_QQ for Rigid Rotor (Option B)")
println("="^70)
println("(Tip: set JULIA_NUM_THREADS to use all cores)")

PARALLEL_SW = get(ENV, "UT_PARALLEL_SW", "0") == "1"
G_SIGN = parse(Float64, get(ENV, "UT_G_SIGN", "1"))

@variables Δ₊ Δ₋ g c₀₁ c₁₂

ω_c = (Δ₊ + Δ₋) / 2
B = (Δ₊ - Δ₋) / 4

L = nlevel_ops(3, :L)

H_rotor = 0 * L[1, 1] + 2B * L[2, 2] + 6B * L[3, 3]
H_cav = ω_c * a'() * a()

cos_theta = c₀₁ * (L[1, 2] + L[2, 1]) + c₁₂ * (L[2, 3] + L[3, 2])
H_int = (G_SIGN * g) * cos_theta * (a'() + a())

H = normal_form(H_rotor + H_cav + H_int)

P_L0 = Subspace(L[1, 1] => 1)

println("Running SW4 with include_QQ=false (legacy)...")
t_false = @elapsed result_false = schrieffer_wolff(
    H,
    P_L0;
    order = 4,
    simplify_mode = :fast,
    include_QQ = false,
    parallel = PARALLEL_SW,
)
println("  done in $(round(t_false, digits=2)) s (parallel=$(PARALLEL_SW))")

println("Running SW4 with include_QQ=true (new)...")
t_true = @elapsed result_true = schrieffer_wolff(
    H,
    P_L0;
    order = 4,
    simplify_mode = :fast,
    include_QQ = true,
    parallel = PARALLEL_SW,
)
println("  done in $(round(t_true, digits=2)) s (parallel=$(PARALLEL_SW))")

println("H_P term counts: false=$(length(result_false.H_P.terms)), true=$(length(result_true.H_P.terms))")

function clean_rationals(expr)
    uexpr = Symbolics.unwrap(expr)
    clean = SymbolicUtils.Postwalk(
        x -> (x isa Rational && isone(denominator(x))) ? numerator(x) : x;
        threaded = false,
    )
    return Symbolics.wrap(clean(uexpr))
end

function extract_hermitian_coeffs(H_eff::QuExpr)
    const_coeff = nothing
    for (term, coeff) in H_eff.terms
        if isempty(term.bares.v)
            const_coeff = coeff
            break
        end
    end

    coeff_a2 = extract_coefficient(H_eff, a()^2; simplify_coeff = false)
    coeff_adag2 = extract_coefficient(H_eff, a'()^2; simplify_coeff = false)
    coeff_adaga = extract_coefficient(H_eff, a'() * a(); simplify_coeff = false)

    coeff_a4 = extract_coefficient(H_eff, a()^4; simplify_coeff = false)
    coeff_adag4 = extract_coefficient(H_eff, a'()^4; simplify_coeff = false)
    coeff_adag3a = extract_coefficient(H_eff, a'()^3 * a(); simplify_coeff = false)
    coeff_adag2a2 = extract_coefficient(H_eff, a'()^2 * a()^2; simplify_coeff = false)

    return (
        E0 = const_coeff,
        A = coeff_adag2,
        Ω = coeff_adaga,
        κ = coeff_adag4,
        μ = coeff_adag3a,
        ν = coeff_adag2a2,
    )
end

println("Extracting coefficients...")
t_coeff_false = @elapsed coeffs_false = extract_hermitian_coeffs(result_false.H_P)
t_coeff_true = @elapsed coeffs_true = extract_hermitian_coeffs(result_true.H_P)
println("  false in $(round(t_coeff_false, digits=2)) s")
println("  true  in $(round(t_coeff_true, digits=2)) s")

function numeric_value(expr, subs)
    expr === nothing && return 0.0
    expr isa Num || return expr
    val = Symbolics.substitute(expr, subs)
    val = val isa Num ? Symbolics.value(val) : val
    return val
end

function near_zero(val; tol = 1e-8)
    if val isa Num
        return false
    end
    if !isfinite(val)
        return false
    end
    return abs(val) < tol
end

function check_diff(name, diff)
    subs_zero_1 = Dict(Δ₊ => 7.0, Δ₋ => 5.0, g => 0.1, c₀₁ => 0.6, c₁₂ => 0.0)
    subs_zero_2 = Dict(Δ₊ => 9.0, Δ₋ => 4.0, g => 0.2, c₀₁ => 0.5, c₁₂ => 0.0)
    subs_nonzero = Dict(Δ₊ => 7.0, Δ₋ => 5.0, g => 0.1, c₀₁ => 0.6, c₁₂ => 0.4)

    v0_1 = numeric_value(diff, subs_zero_1)
    v0_2 = numeric_value(diff, subs_zero_2)
    v1 = numeric_value(diff, subs_nonzero)

    zero_ok = near_zero(v0_1) && near_zero(v0_2)
    nonzero_ok = !near_zero(v1)

    println("  $name: c12=0 -> $(zero_ok ? "0" : "NONZERO"), c12=0.4 -> $(nonzero_ok ? "nonzero" : "0")")
    if !zero_ok
        println("    values: c12=0 -> ($v0_1, $v0_2)")
    end
end

println("\nComparing coefficients (true - false):")
diffs = (
    E0 = coeffs_true.E0 - coeffs_false.E0,
    A = coeffs_true.A - coeffs_false.A,
    Ω = coeffs_true.Ω - coeffs_false.Ω,
    κ = coeffs_true.κ - coeffs_false.κ,
    μ = coeffs_true.μ - coeffs_false.μ,
    ν = coeffs_true.ν - coeffs_false.ν,
)

for (name, diff) in pairs(diffs)
    check_diff(string(name), diff)
end

println("\nIf any coefficient is NONZERO at c12=0, include_QQ is changing c01-only terms (bug).")
