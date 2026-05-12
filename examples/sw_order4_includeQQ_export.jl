using UnitaryTransformations
using QuantumAlgebra
using Symbolics

L = nlevel_ops(3, :L)
@variables Δ δ ω g1 g2

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

sections = [
    "Input Hamiltonian" => "H = " * string(H),
    "Generator" => "S = " * string(result.S),
    "Effective Hamiltonian" => "H_{\\mathrm{eff}} = " * string(result.H_eff),
    "Projected Hamiltonian" => "H_P = " * string(result.H_P),
]

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

tex_path = joinpath(outdir, "sw_order4_includeQQ_full.tex")
write_expression_dump(
    tex_path,
    sections;
    mode = :latex,
    tex = true,
    header = [
        "% 3-level ladder + cavity SW to 4th order (include_QQ=true)",
        "% Variables: \\Delta, \\delta, \\omega, g_1, g_2",
    ],
)

unicode_path = joinpath(outdir, "sw_order4_includeQQ_full.unicode.txt")
write_expression_dump(
    unicode_path,
    sections;
    mode = :unicode,
    tex = false,
    header = [
        "# 3-level ladder + cavity SW to 4th order (include_QQ=true)",
        "# Unicode mode: keeps Greek letters and uses a'()",
    ],
)

println("Wrote ", tex_path)
println("Wrote ", unicode_path)
