#=
Compare SW cavity H_eff vs ED ground energy
===========================================

Uses the full SW cavity Hamiltonian (Option B, L=0) and compares its ground
state energy to ED of the full rotor+cavity model. Plots SW2 vs SW4 errors.
=#

using Serialization
using LinearAlgebra
using SparseArrays

using UnitaryTransformations
using UnitaryTransformations: extract_coefficient
using QuantumAlgebra
using Symbolics
using KrylovKit
using PyPlot

include("rigid_rotor_ed.jl")

const CACHE_DIR = joinpath(@__DIR__, "cache")
const PLOT_DIR = joinpath(@__DIR__, "plots")
const PLOT_FILE = joinpath(PLOT_DIR, "ed_vs_sw_ground.png")

const SW_BASIS = get(ENV, "UT_SW_BASIS", "detuning")
const SW4_DIAGONAL_ONLY = get(ENV, "UT_SW4_DIAGONAL_ONLY", SW_BASIS == "detuning" ? "0" : "1") == "1"
const INCLUDE_QQ = get(ENV, "UT_INCLUDE_QQ", "1") == "1"
const USE_FIXED_C = get(ENV, "UT_USE_FIXED_C", "1") == "1"
const G_SIGN = parse(Float64, get(ENV, "UT_G_SIGN", "1"))
const SHOW_PLOT = get(ENV, "UT_SHOW_PLOT", "1") == "1"

const CACHE_FILE = joinpath(
    CACHE_DIR,
    "sw_cavity_expr_$(SW_BASIS)_diag$(SW4_DIAGONAL_ONLY)_qq$(INCLUDE_QQ)_gsign$(G_SIGN).jls",
)

mkpath(CACHE_DIR)
mkpath(PLOT_DIR)

if SHOW_PLOT
    PyPlot.ion()
end

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

function compute_sw_cavity_expr(
    ;
    include_QQ::Bool,
    sw4_diagonal_only::Bool,
    basis::String,
    g_sign::Float64,
)
    L = nlevel_ops(3, :L)
    if basis == "detuning"
        @variables Dp Dm g c01 c12
        wc = (Dp + Dm) / 2
        B = (Dp - Dm) / 4
        vars = (Dp = Dp, Dm = Dm, g = g, c01 = c01, c12 = c12)
    else
        @variables B wc g c01 c12
        vars = (B = B, wc = wc, g = g, c01 = c01, c12 = c12)
    end

    H_rotor = 0 * L[1, 1] + 2B * L[2, 2] + 6B * L[3, 3]
    H_cav = wc * a'() * a()
    cos_theta = c01 * (L[1, 2] + L[2, 1]) + c12 * (L[2, 3] + L[3, 2])
    H_int = (g_sign * g) * cos_theta * (a'() + a())
    H = normal_form(H_rotor + H_cav + H_int)

    P_L0 = Subspace(L[1, 1] => 1)

    result_o2 = schrieffer_wolff(
        H,
        P_L0;
        order = 2,
        simplify_mode = :fractions,
        include_QQ = include_QQ,
        parallel = false,
    )

    sw4_simplify_mode = :fast
    if !sw4_diagonal_only && basis == "detuning"
        sw4_simplify_mode = :none
    end

    result_o4 = schrieffer_wolff(
        H,
        P_L0;
        order = 4,
        simplify_mode = sw4_simplify_mode,
        include_QQ = include_QQ,
        diagonal_only = sw4_diagonal_only,
        parallel = false,
    )

    return (
        include_QQ = include_QQ,
        sw4_diagonal_only = sw4_diagonal_only,
        basis = basis,
        coeffs_o2 = extract_cavity_coeffs(result_o2.H_P),
        coeffs_o4 = extract_cavity_coeffs(result_o4.H_P),
        vars = vars,
    )
end

function load_sw_cache(
    cache_file;
    include_QQ::Bool,
    sw4_diagonal_only::Bool,
    basis::String,
    g_sign::Float64,
)
    if isfile(cache_file)
        cache = deserialize(cache_file)
        if hasproperty(cache, :include_QQ) &&
           hasproperty(cache, :sw4_diagonal_only) &&
           hasproperty(cache, :basis) &&
           cache.include_QQ == include_QQ &&
           cache.sw4_diagonal_only == sw4_diagonal_only &&
           cache.basis == basis
            return cache, true
        end
    end
    cache = compute_sw_cavity_expr(
        ;
        include_QQ = include_QQ,
        sw4_diagonal_only = sw4_diagonal_only,
        basis = basis,
        g_sign = g_sign,
    )
    open(cache_file, "w") do io
        serialize(io, cache)
    end
    return cache, false
end

function eval_expr(expr, subs)
    expr === nothing && return 0.0
    val = Symbolics.substitute(expr, subs)
    val = val isa Num ? Symbolics.value(val) : val
    return Float64(val)
end

function eval_coeffs(coeffs::NamedTuple, subs)
    return (
        E0 = eval_expr(coeffs.E0, subs),
        A = eval_expr(coeffs.A, subs),
        Ω = eval_expr(coeffs.Ω, subs),
        κ = eval_expr(coeffs.κ, subs),
        μ = eval_expr(coeffs.μ, subs),
        ν = eval_expr(coeffs.ν, subs),
    )
end

function build_ed_operators(
    omega_c::Float64,
    B::Float64;
    l_max::Int,
    n_max::Int,
    gauge_term::Bool,
    cosθ::SparseMatrixCSC{Float64,Int},
)
    cosθ² = cosθ * cosθ

    l_all = 0:l_max
    L2 = spdiagm(0 => l_all .* (l_all .+ 1))
    a = gen_a(n_max)

    id_ph = spdiagm(0 => ones(n_max + 1))
    id_l = spdiagm(0 => ones(l_max + 1))

    H0 = B * kron(id_ph, L2) + omega_c * kron(a' * a, id_l)
    V = kron(a + a', cosθ)
    W = gauge_term ? kron(id_ph, cosθ²) : spzeros(size(H0))

    return H0, V, W
end

function build_cavity_matrix(coeffs::NamedTuple, n_max::Int)
    a = gen_a(n_max)
    adag = a'
    id_ph = spdiagm(0 => ones(n_max + 1))

    H = coeffs.E0 * id_ph
    H += coeffs.A * (adag^2 + a^2)
    H += coeffs.Ω * (adag * a)
    H += coeffs.κ * (adag^4 + a^4)
    H += coeffs.μ * (adag^3 * a + adag * a^3)
    H += coeffs.ν * (adag^2 * a^2)

    return H
end

function ground_state_ed(H)
    vals, vecs, _ = eigsolve(H, 1, :SR; ishermitian = true)
    return real(vals[1]), vecs[1]
end

function ground_energy_ed(H)
    val, _ = ground_state_ed(H)
    return val
end

println("="^70)
println("  Compare SW cavity H_eff vs ED ground energy")
println("="^70)

# ---- Parameters ----
B = 1.0
omega_c = parse(Float64, get(ENV, "UT_OMEGA_C", "0.2"))
g_max = parse(Float64, get(ENV, "UT_G_MAX", "1.0"))
num_points = parse(Int, get(ENV, "UT_NUM_POINTS", "51"))
l_max = parse(Int, get(ENV, "UT_L_MAX", "2"))
n_max = parse(Int, get(ENV, "UT_N_MAX", "80"))
n_eff = parse(Int, get(ENV, "UT_N_EFF", string(n_max)))
gauge_term = get(ENV, "UT_GAUGE_TERM", "0") == "1"

if USE_FIXED_C
    c01 = 1 / sqrt(3)
    c12 = sqrt(2 / 15)
else
    cosθ_ref = gen_cosθ(2)
    c01 = cosθ_ref[1, 2]
    c12 = cosθ_ref[2, 3]
end

println("Using c01=$(c01), c12=$(c12)")
println("omega_c=$(omega_c), B=$(B)")
println("g range: 0 → $(g_max) with $(num_points) points")
println("n_max=$(n_max), n_eff=$(n_eff), l_max=$(l_max)")

g_values = range(0.0, g_max; length = num_points)

cache, loaded = load_sw_cache(
    CACHE_FILE;
    include_QQ = INCLUDE_QQ,
    sw4_diagonal_only = SW4_DIAGONAL_ONLY,
    basis = SW_BASIS,
    g_sign = G_SIGN,
)
println(loaded ? "Loaded cached SW expressions." : "Computed and cached SW expressions.")
println("SW basis=$(SW_BASIS), include_QQ=$(INCLUDE_QQ), SW4 diagonal_only=$(SW4_DIAGONAL_ONLY)")
println("Coupling sign: $(G_SIGN)")

subs_base = Dict(
    cache.vars.c01 => c01,
    cache.vars.c12 => c12,
)
if cache.basis == "detuning"
    Dp = omega_c + 2B
    Dm = omega_c - 2B
    subs_base[cache.vars.Dp] = Dp
    subs_base[cache.vars.Dm] = Dm
else
    subs_base[cache.vars.B] = B
    subs_base[cache.vars.wc] = omega_c
end

sw2 = zeros(length(g_values))
sw4 = zeros(length(g_values))

for (i, g_val) in enumerate(g_values)
    subs = copy(subs_base)
    subs[cache.vars.g] = g_val
    coeffs_o2 = eval_coeffs(cache.coeffs_o2, subs)
    coeffs_o4 = eval_coeffs(cache.coeffs_o4, subs)
    H2 = build_cavity_matrix(coeffs_o2, n_eff)
    H4 = build_cavity_matrix(coeffs_o4, n_eff)
    sw2[i] = ground_energy_ed(H2)
    sw4[i] = ground_energy_ed(H4)
end

println("Building ED operators...")
cosθ_ed = gen_cosθ(l_max)
if l_max == 2
    cosθ_ed = spzeros(Float64, l_max + 1, l_max + 1)
    cosθ_ed[1, 2] = c01
    cosθ_ed[2, 1] = c01
    cosθ_ed[2, 3] = c12
    cosθ_ed[3, 2] = c12
end

H0, V, W = build_ed_operators(
    omega_c,
    B;
    l_max = l_max,
    n_max = n_max,
    gauge_term = gauge_term,
    cosθ = cosθ_ed,
)

ed = zeros(length(g_values))
for (i, g_val) in enumerate(g_values)
    g_eff = G_SIGN * g_val
    H = H0 + g_eff * V + (g_eff^2 / omega_c) * W
    ed[i] = ground_energy_ed(H)
end

# Check photon cutoff at largest coupling
g_max_eff = G_SIGN * g_max
H_max = H0 + g_max_eff * V + (g_max_eff^2 / omega_c) * W
E_max, psi_max = ground_state_ed(H_max)
psi_mat = reshape(psi_max, l_max + 1, n_max + 1)
weights_n = vec(sum(abs2, psi_mat; dims = 1))
boundary_weight = weights_n[end]
mean_n = sum((0:n_max) .* weights_n)
println("Photon cutoff check at g_max: <n>=$(round(mean_n, digits=4)), weight(n_max)=$(boundary_weight)")
if boundary_weight > 1e-6
    println("  ⚠ Consider increasing n_max (boundary weight > 1e-6)")
end

err2 = abs.(sw2 .- ed)
err4 = abs.(sw4 .- ed)
better_count = count(err4 .< err2)

println("Max error SW2: ", maximum(err2))
println("Max error SW4: ", maximum(err4))
println("SW4 better at $(better_count) / $(length(err2)) points")

fig, axes = subplots(2, 1; figsize = (6, 7), sharex = true)

axes[1].plot(g_values, ed; label = "ED", color = "black", linewidth = 2)
axes[1].plot(g_values, sw2; label = "SW2", linestyle = "--")
axes[1].plot(g_values, sw4; label = "SW4", linestyle = "-")
axes[1].set_ylabel("Ground energy")
axes[1].legend()
axes[1].set_ylim(-0.5,0.001)

axes[2].plot(g_values, err2; label = "|SW2 - ED|", linestyle = "--")
axes[2].plot(g_values, err4; label = "|SW4 - ED|", linestyle = "-")
axes[2].set_yscale("log")
axes[2].set_xlabel("g")
axes[2].set_ylabel("Abs error")
axes[2].legend()

fig.tight_layout()

savefig(PLOT_FILE; dpi = 200)


println("Saved plot to $(PLOT_FILE)")
if SHOW_PLOT
    plt.show()
end
