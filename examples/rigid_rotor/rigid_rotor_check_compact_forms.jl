#=
Check compact x=Δ-/Δ+ forms for κ, μ, ν
========================================

Validates the compact rational forms against the SW coefficients obtained
from the full detuning-basis expansion (Option B, L=0).
=#

using Serialization
using Random

using UnitaryTransformations
using UnitaryTransformations: extract_coefficient
using QuantumAlgebra
using Symbolics

const CACHE_FILE = joinpath(@__DIR__, "cache", "sw_cavity_expr_detuning_diag0_qq1_gsign1.jls")

function extract_constant(expr::QuExpr)
    for (term, coeff) in expr.terms
        if isempty(term.bares.v)
            full_coeff = coeff
            for p in term.params
                full_coeff = full_coeff * UnitaryTransformations.param_to_symbolic(p)
            end
            return full_coeff
        end
    end
    return 0
end

function extract_cavity_coeffs(H_eff::QuExpr)
    return (
        E0 = extract_constant(H_eff),
        A = extract_coefficient(H_eff, a'()^2; simplify_coeff = false),
        Ω = extract_coefficient(H_eff, a'() * a(); simplify_coeff = false),
        κ = extract_coefficient(H_eff, a'()^4; simplify_coeff = false),
        μ = extract_coefficient(H_eff, a'()^3 * a(); simplify_coeff = false),
        ν = extract_coefficient(H_eff, a'()^2 * a()^2; simplify_coeff = false),
    )
end

function compute_sw_coeffs()
    L = nlevel_ops(3, :L)
    @variables Dp Dm g c01 c12
    wc = (Dp + Dm) / 2
    B = (Dp - Dm) / 4

    H_rotor = 0 * L[1, 1] + 2B * L[2, 2] + 6B * L[3, 3]
    H_cav = wc * a'() * a()
    cos_theta = c01 * (L[1, 2] + L[2, 1]) + c12 * (L[2, 3] + L[3, 2])
    H_int = g * cos_theta * (a'() + a())
    H = normal_form(H_rotor + H_cav + H_int)

    P_L0 = Subspace(L[1, 1] => 1)
    result_o4 = schrieffer_wolff(
        H,
        P_L0;
        order = 4,
        simplify_mode = :fast,
        include_QQ = true,
        diagonal_only = false,
        parallel = false,
    )

    return (
        coeffs_o4 = extract_cavity_coeffs(result_o4.H_P),
        vars = (Dp = Dp, Dm = Dm, g = g, c01 = c01, c12 = c12),
    )
end

function load_sw_coeffs()
    if isfile(CACHE_FILE)
        return deserialize(CACHE_FILE)
    end
    cache = compute_sw_coeffs()
    open(CACHE_FILE, "w") do io
        serialize(io, cache)
    end
    return cache
end

function eval_expr(expr, subs)
    expr === nothing && return 0.0
    val = Symbolics.substitute(expr, subs)
    val = val isa Num ? Symbolics.value(val) : val
    return BigFloat(val)
end

function poly_eval(coeffs, x)
    result = zero(x)
    power = one(x)
    for c in coeffs
        result += c * power
        power *= x
    end
    return result
end

function pal_poly(coeffs, x)
    full = vcat(coeffs, reverse(coeffs[1:end-1]))
    return poly_eval(full, x)
end

function kappa_formula(Dp, Dm, g, c01, c12)
    x = Dm / Dp
    pκ = BigFloat.([
        -390,
        10937,
        -122714,
        674562,
        -1612057,
        -491071,
        9461936,
        -10654341,
        -8219134,
        4444439,
        13435616,
        -6873027,
        -1205945,
        1505781,
        -397949,
        45312,
        -1953,
    ])
    qκ = BigFloat.([
        0,
        -781,
        17890,
        -153921,
        558762,
        -339119,
        -2748101,
        4411233,
        2771667,
        -1335087,
        -5227870,
        1634703,
        777812,
        -433984,
        70703,
        -3906,
    ])
    numer = c01^4 * poly_eval(pκ, x) + c01^2 * c12^2 * poly_eval(qκ, x)
    denom = x^2 * (5 - x)^5 * (x + 2)^2 * (2x + 1)^2 * (5x - 1)^6
    return (16 * g^4 / Dp^3) * numer / denom
end

function mu_formula(Dp, Dm, g, c01, c12)
    x = Dm / Dp
    pμ = BigFloat.([
        -81091461181,
        4784396209717,
        -130155849769592,
        2161173397439575,
        -24415642023368225,
        197709668281328466,
        -1173546944400198780,
        5090610808136502772,
        -15428354239729496356,
        26945226505785849988,
        9740945307635679417,
        -235255289959426826036,
        825831847158476295675,
        -1655421335785389783999,
        1760776773890129480506,
        707358223608740601279,
        -7271370234578271913529,
        17125649636140287001396,
        -26396872845443003702588,
        30227906027567840262997,
    ])
    qμ = BigFloat.([
        -216243896484,
        11623109436035,
        -285825776275634,
        4250185125402832,
        -42499926683349609,
        299891522238296484,
        -1515379300946658984,
        5364724543521525703,
        -11895278791979861437,
        6907134017956398188,
        60262302777292632178,
        -268108549639864182271,
        603861902214556427069,
        -744477593756820620657,
        -137151241914754094,
        2246841393679613753290,
        -5756512261690495069553,
        9116008724734496141949,
        -10513722669814370690000,
    ])
    Pμ = pal_poly(pμ, x)
    Qμ = x * pal_poly(qμ, x)
    numer = c01^4 * Pμ + c01^2 * c12^2 * Qμ
    denom = x^3 * (1 - x)^11 * (5 - x)^10 * (x + 2)^2 * (2x + 1)^2 * (5x - 1)^10
    return (128 * g^4 / (531441 * Dp^3)) * numer / denom
end

function nu_formula(Dp, Dm, g, c01, c12)
    x = Dm / Dp
    pν = BigFloat.([
        5911567520141601,
        -481201596139526367,
        18646502734731445312,
        -458005314578593139648,
        8011712841756173217773,
        -106324076148778660341797,
        1113730129873134503999999,
        -9460065916572125518585547,
        66443295911081766520160742,
        -391606102605085161921316054,
        1959214180946961654149255663,
        -8397484281067667957041679287,
        31068849798981040111207124098,
        -99846450574215727490909924907,
        280195742376855922921192632529,
        -689682506449779255161141741199,
        1494633437911217399009252925148,
        -2860851126627629007166893255149,
        4849206878264122173844812779794,
        -7294226447811573907507414076289,
        9752759009067135315798844690865,
        -11604241845603451351953818999217,
        12295512675286272430358288423906,
    ])
    qν = BigFloat.([
        15764180053710937,
        -1196501266076660156,
        43106204506069335937,
        -981466499934404296875,
        15866343795122384765625,
        -194001069241740662929687,
        1866610686256322648765625,
        -14520189471659502499321875,
        93127414515052683970617187,
        -499834946911799397645817031,
        2271365489045527262567829787,
        -8821586351582792570892355740,
        29510370284824999048032484524,
        -85583660762548493701664849082,
        216360389893368430771174822572,
        -479019947554053131089732783884,
        932437610394442360474829539275,
        -1600997220704536567752397514044,
        2431145318235870115311220603611,
        -3271645610462915394861297635226,
        3907395613867381907476398386742,
        -4145266204287696662684153396445,
    ])
    Pν = pal_poly(pν, x)
    Qν = x * pal_poly(qν, x)
    numer = c01^4 * Pν + c01^2 * c12^2 * Qν
    denom = x^3 * (5 - x)^12 * (1 - x)^17 * (5x - 1)^12
    return -(16 * g^4 / (387420489 * Dp^3)) * numer / denom
end

function rel_err(a, b)
    return abs(a - b) / max(BigFloat(1), abs(a))
end

println("="^70)
println("  Check compact x=Δ-/Δ+ forms (κ, μ, ν)")
println("="^70)

setprecision(256) do
    cache = load_sw_coeffs()
    coeffs = cache.coeffs_o4
    vars = cache.vars

    c01 = BigFloat(1 / sqrt(3))
    c12 = BigFloat(sqrt(2 / 15))

    samples = [
        (Dp = BigFloat(2.2), Dm = BigFloat(-1.8), g = BigFloat(0.08)),
        (Dp = BigFloat(3.1), Dm = BigFloat(1.7), g = BigFloat(0.05)),
        (Dp = BigFloat(5.0), Dm = BigFloat(2.0), g = BigFloat(0.1)),
        (Dp = BigFloat(2.5), Dm = BigFloat(-0.7), g = BigFloat(0.12)),
    ]

    Random.seed!(42)
    for _ in 1:2
        Dp = BigFloat(2.0 + rand() * 3.0)
        Dm = BigFloat(-1.5 + rand() * 2.5)
        g = BigFloat(0.03 + rand() * 0.07)
        push!(samples, (Dp = Dp, Dm = Dm, g = g))
    end

    for (i, s) in enumerate(samples)
        subs = Dict(
            vars.Dp => s.Dp,
            vars.Dm => s.Dm,
            vars.g => s.g,
            vars.c01 => c01,
            vars.c12 => c12,
        )

        sw_k = eval_expr(coeffs.κ, subs)
        sw_m = eval_expr(coeffs.μ, subs)
        sw_n = eval_expr(coeffs.ν, subs)

        f_k = kappa_formula(s.Dp, s.Dm, s.g, c01, c12)
        f_m = mu_formula(s.Dp, s.Dm, s.g, c01, c12)
        f_n = nu_formula(s.Dp, s.Dm, s.g, c01, c12)

        println("Sample $(i): Dp=$(s.Dp), Dm=$(s.Dm), g=$(s.g)")
        println("  κ: rel_err=", rel_err(sw_k, f_k))
        println("  μ: rel_err=", rel_err(sw_m, f_m))
        println("  ν: rel_err=", rel_err(sw_n, f_n))
    end
end
