using UnitaryTransformations
using QuantumAlgebra
using QuantumAlgebra: QuTerm
using Symbolics
using Symbolics: substitute, simplify

function map_expr_coefficients(expr::QuExpr, f)
    out = Dict{QuTerm,Number}()
    for (term, coeff) in expr.terms
        out[term] = f(coeff)
    end
    return normal_form(QuExpr(out))
end

L = nlevel_ops(3, :L)
@variables Δ δ ω g1 g2 d1 d3

H0 = 0 * L[1, 1] + Δ * L[2, 2] + (Δ + δ) * L[3, 3] + ω * a'() * a()
V =
    g1 * (L[1, 2] * a'() + L[2, 1] * a()) +
    g2 * (L[2, 3] * a'() + L[3, 2] * a())
H = normal_form(H0 + V)

P = Subspace(L[1, 1] => 1)
result = schrieffer_wolff(
    H,
    P;
    order = 4,
    include_QQ = true,
    diagonal_only = false,
    simplify_mode = :none,
)

detuning_subs = Dict(
    Δ => d1 + ω,
    δ => d3 - d1 + ω,
)

function sw_scalar_simplify(c)
    c1 = substitute(c, detuning_subs)
    c2 = simplify(c1; expand = true)
    c3 = simplify(c2; expand = false)
    return c3
end

H_eff_detuned = map_expr_coefficients(result.H_eff, sw_scalar_simplify)
H_P_detuned = map_expr_coefficients(result.H_P, sw_scalar_simplify)

sections_tex = [
    "Input Hamiltonian" => "H = " * string(H),
    "Generator" => "S = " * string(result.S),
    "Effective Hamiltonian (raw)" => "H_{\\mathrm{eff}} = " * string(result.H_eff),
    "Effective Hamiltonian (detuning vars)" => "H_{\\mathrm{eff}} = " * string(H_eff_detuned),
    "Projected Hamiltonian (raw)" => "H_P = " * string(result.H_P),
    "Projected Hamiltonian (detuning vars)" => "H_P = " * string(H_P_detuned),
]

sections_unicode = [
    "Input Hamiltonian" => "H = " * string(H),
    "Generator" => "S = " * string(result.S),
    "Effective Hamiltonian (raw)" => "H_eff = " * string(result.H_eff),
    "Effective Hamiltonian (detuning vars)" => "H_eff = " * string(H_eff_detuned),
    "Projected Hamiltonian (raw)" => "H_P = " * string(result.H_P),
    "Projected Hamiltonian (detuning vars)" => "H_P = " * string(H_P_detuned),
]

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

tex_path = joinpath(outdir, "sw_order4_includeQQ_full.tex")
write_expression_dump(
    tex_path,
    sections_tex;
    mode = :latex,
    tex = true,
    header = [
        "% 3-level ladder + cavity SW to 4th order (include_QQ=true)",
        "% Variables: \\Delta, \\delta, \\omega, g_1, g_2",
        "% Detunings: d_1 = \\Delta - \\omega, d_3 = \\Delta + \\delta - 2\\omega",
        "% Substitution: \\Delta => d_1 + \\omega, \\delta => d_3 - d_1 + \\omega",
    ],
)

unicode_path = joinpath(outdir, "sw_order4_includeQQ_full.unicode.txt")
write_expression_dump(
    unicode_path,
    sections_unicode;
    mode = :unicode,
    tex = false,
    header = [
        "# 3-level ladder + cavity SW to 4th order (include_QQ=true)",
        "# Unicode mode: keeps Greek letters and uses a'",
        "# Detunings: d1 = Δ - ω, d3 = Δ + δ - 2ω",
        "# Substitution: Δ => d1 + ω, δ => d3 - d1 + ω",
    ],
)

hp_detuned_str = string(H_P_detuned)
println("Detuned H_P contains d1: ", occursin("d1", hp_detuned_str))
println("Detuned H_P contains d3: ", occursin("d3", hp_detuned_str))
println(
    "Detuned H_P contains mixed g1^2 g2^2 structure: ",
    occursin("g1^2", hp_detuned_str) && occursin("g2^2", hp_detuned_str),
)

println("Wrote ", tex_path)
println("Wrote ", unicode_path)
