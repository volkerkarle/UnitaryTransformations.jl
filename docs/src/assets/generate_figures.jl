#=
Generate Publication-Ready Figures for UnitaryTransformations.jl Documentation

This script creates figures comparing ground state energies at different
Schrieffer-Wolff orders as a function of coupling strength g.

Figures generated:
1. Two-level system: E_g vs g/Δ comparing SW orders 2,4 with exact solution
2. Jaynes-Cummings: Ground state energy shift vs g/Δ at orders 2,4
3. Error analysis: Relative error of SW approximations vs coupling strength
=#

using CairoMakie
using LaTeXStrings

# Set up publication-quality theme
function set_publication_theme!()
    theme = Theme(
        # Font settings
        fontsize = 14,
        font = "CMU Serif",
        
        # Axis settings
        Axis = (
            xlabelsize = 16,
            ylabelsize = 16,
            xticklabelsize = 12,
            yticklabelsize = 12,
            titlesize = 16,
            xgridvisible = true,
            ygridvisible = true,
            xgridstyle = :dash,
            ygridstyle = :dash,
            xgridcolor = (:gray, 0.3),
            ygridcolor = (:gray, 0.3),
            spinewidth = 1.5,
            xtickwidth = 1.5,
            ytickwidth = 1.5,
            xminorticksvisible = true,
            yminorticksvisible = true,
        ),
        
        # Legend settings  
        Legend = (
            framevisible = true,
            framewidth = 1,
            backgroundcolor = :white,
            labelsize = 12,
            padding = (8, 8, 6, 6),
        ),
        
        # Line settings
        Lines = (
            linewidth = 2.5,
        ),
        
        # Scatter settings
        Scatter = (
            markersize = 8,
        ),
        
        # Figure settings
        figure_padding = 10,
        backgroundcolor = :white,
    )
    set_theme!(theme)
end

# Color palette (colorblind-friendly)
const COLORS = [
    colorant"#0077BB",  # Blue
    colorant"#EE7733",  # Orange  
    colorant"#009988",  # Teal
    colorant"#CC3311",  # Red
    colorant"#EE3377",  # Magenta
    colorant"#33BBEE",  # Cyan
]

"""
    figure1_two_level_system(; save_path=nothing)

Create figure comparing ground state energy of a two-level system
H = Δ/2 σz + g σx at different SW orders vs exact solution.

Exact: E_g = -√(Δ²/4 + g²)
SW order 2: E_g ≈ -Δ/2 - g²/Δ  
SW order 4: E_g ≈ -Δ/2 - g²/Δ + g⁴/Δ³
"""
function figure1_two_level_system(; save_path=nothing)
    Δ = 1.0  # Energy scale
    g_range = range(0, 0.8, length=100)
    
    # Compute energies using the actual SW results (evaluated symbolically):
    # Order 2: H_P = (Δ² + 2g²) / (-2Δ) = -Δ/2 - g²/Δ
    # Order 4: H_P = (Δ⁴/2 + Δ²g² - g⁴) / (-Δ³) = -Δ/2 - g²/Δ + g⁴/Δ³
    # Order 6: H_P = (Δ⁶/2 + Δ⁴g² - Δ²g⁴ + 2g⁶) / (-Δ⁵) = -Δ/2 - g²/Δ + g⁴/Δ³ - 2g⁶/Δ⁵
    E_exact = [-sqrt(Δ^2/4 + g^2) for g in g_range]
    E_SW2 = [-Δ/2 - g^2/Δ for g in g_range]
    E_SW4 = [-Δ/2 - g^2/Δ + g^4/Δ^3 for g in g_range]
    E_SW6 = [-Δ/2 - g^2/Δ + g^4/Δ^3 - 2*g^6/Δ^5 for g in g_range]
    
    fig = Figure(size = (700, 500))
    
    # Main plot
    ax1 = Axis(fig[1, 1],
        xlabel = L"g / \Delta",
        ylabel = L"E_g / \Delta",
        title = L"\mathrm{Two-level system: } H = \frac{\Delta}{2}\sigma_z + g \sigma_x",
    )
    
    # Plot lines
    lines!(ax1, g_range, E_exact, color = :black, linewidth = 3, 
           linestyle = :solid, label = "Exact")
    lines!(ax1, g_range, E_SW2, color = COLORS[1], linewidth = 2.5,
           linestyle = :dash, label = "SW order 2")
    lines!(ax1, g_range, E_SW4, color = COLORS[2], linewidth = 2.5,
           linestyle = :dashdot, label = "SW order 4")
    lines!(ax1, g_range, E_SW6, color = COLORS[3], linewidth = 2.5,
           linestyle = :dot, label = "SW order 6")
    
    # Add legend
    axislegend(ax1, position = :lb)
    
    # Inset: relative error
    ax2 = Axis(fig[1, 1],
        width = Relative(0.4),
        height = Relative(0.35),
        halign = 0.95,
        valign = 0.95,
        xlabel = L"g / \Delta",
        ylabel = "Relative error (%)",
        xlabelsize = 11,
        ylabelsize = 11,
        xticklabelsize = 9,
        yticklabelsize = 9,
        backgroundcolor = (:white, 0.9),
    )
    
    # Relative errors (avoid division by very small numbers)
    mask = abs.(E_exact) .> 0.01
    err_SW2 = 100 .* abs.(E_SW2[mask] .- E_exact[mask]) ./ abs.(E_exact[mask])
    err_SW4 = 100 .* abs.(E_SW4[mask] .- E_exact[mask]) ./ abs.(E_exact[mask])
    err_SW6 = 100 .* abs.(E_SW6[mask] .- E_exact[mask]) ./ abs.(E_exact[mask])
    
    lines!(ax2, g_range[mask], err_SW2, color = COLORS[1], linewidth = 2)
    lines!(ax2, g_range[mask], err_SW4, color = COLORS[2], linewidth = 2)
    lines!(ax2, g_range[mask], err_SW6, color = COLORS[3], linewidth = 2)
    
    ylims!(ax2, 0, 15)
    
    if save_path !== nothing
        save(save_path, fig, px_per_unit = 3)
        println("Saved: $save_path")
    end
    
    return fig
end

"""
    figure2_jaynes_cummings(; save_path=nothing)

Create figure showing the dispersive shift and Kerr nonlinearity
in the Jaynes-Cummings model at different SW orders.

H = ω_c a†a + Δ/2 σz + g(a†σ⁻ + a σ⁺)

Shows:
- Panel (a): Effective cavity frequency shift vs g/Δ
- Panel (b): Kerr coefficient K vs g/Δ (only appears at order 4+)
"""
function figure2_jaynes_cummings(; save_path=nothing)
    Δ = 1.0
    g_range = range(0, 0.5, length=100)
    
    # Dispersive shift: χ = -g²/Δ + higher order corrections
    # At order 2: χ₂ = -g²/Δ
    # At order 4: χ₄ = -g²/Δ + O(g⁴) corrections
    χ_SW2 = [-g^2/Δ for g in g_range]
    # From the SW result, the order-4 correction modifies the a†a coefficient
    χ_SW4 = [-g^2/Δ + 4*g^4/(3*Δ^3) for g in g_range]
    
    # Kerr coefficient K (coefficient of a†²a²)
    # Only appears at order 4: K = 4g⁴/(3Δ³)
    K_SW4 = [4*g^4/(3*Δ^3) for g in g_range]
    
    fig = Figure(size = (900, 400))
    
    # Panel (a): Dispersive shift
    ax1 = Axis(fig[1, 1],
        xlabel = L"g / \Delta",
        ylabel = L"\chi / \Delta",
        title = "(a) Dispersive shift",
    )
    
    lines!(ax1, g_range, χ_SW2, color = COLORS[1], linewidth = 2.5,
           linestyle = :dash, label = L"\mathrm{SW order 2: } -g^2/\Delta")
    lines!(ax1, g_range, χ_SW4, color = COLORS[2], linewidth = 2.5,
           linestyle = :solid, label = "SW order 4")
    
    # Reference: perturbative limit
    lines!(ax1, g_range, -g_range.^2, color = :gray, linewidth = 1.5,
           linestyle = :dot, label = L"\mathrm{Leading order: } -g^2")
    
    axislegend(ax1, position = :rt)
    
    # Panel (b): Kerr coefficient
    ax2 = Axis(fig[1, 2],
        xlabel = L"g / \Delta",
        ylabel = L"K / \Delta",
        title = L"\mathrm{(b) Kerr nonlinearity } K (a^\dagger a)^2",
    )
    
    lines!(ax2, g_range, zeros(length(g_range)), color = COLORS[1], linewidth = 2.5,
           linestyle = :dash, label = L"\mathrm{SW order 2: } K = 0")
    lines!(ax2, g_range, K_SW4, color = COLORS[2], linewidth = 2.5,
           linestyle = :solid, label = L"\mathrm{SW order 4: } K = \frac{4g^4}{3\Delta^3}")
    
    axislegend(ax2, position = :lt)
    
    if save_path !== nothing
        save(save_path, fig, px_per_unit = 3)
        println("Saved: $save_path")
    end
    
    return fig
end

"""
    figure3_convergence(; save_path=nothing)

Create figure showing convergence of SW expansion with increasing order.
Demonstrates the asymptotic nature of the perturbation series.
"""
function figure3_convergence(; save_path=nothing)
    Δ = 1.0
    
    # Different coupling regimes
    g_values = [0.1, 0.3, 0.5, 0.7]
    orders = 2:2:10
    
    fig = Figure(size = (800, 500))
    
    ax = Axis(fig[1, 1],
        xlabel = "SW order n",
        ylabel = L"|E_g^{(n)} - E_g^\mathrm{exact}| / |E_g^\mathrm{exact}|",
        title = "Convergence of SW expansion for two-level system",
        yscale = log10,
        xticks = collect(orders),
    )
    
    for (i, g) in enumerate(g_values)
        E_exact = -sqrt(Δ^2/4 + g^2)
        
        errors = Float64[]
        for n in orders
            # Compute SW approximation to order n
            E_SW = compute_two_level_energy(Δ, g, n)
            push!(errors, abs(E_SW - E_exact) / abs(E_exact))
        end
        
        scatterlines!(ax, collect(orders), errors, 
                     color = COLORS[i], linewidth = 2, markersize = 10,
                     label = L"g/\Delta = %$(g)")
    end
    
    # Add horizontal line for machine precision
    hlines!(ax, [1e-15], color = :gray, linestyle = :dash, linewidth = 1)
    text!(ax, 9, 3e-15, text = "Machine precision", fontsize = 10, color = :gray)
    
    axislegend(ax, position = :rt)
    
    ylims!(ax, 1e-16, 1)
    
    if save_path !== nothing
        save(save_path, fig, px_per_unit = 3)
        println("Saved: $save_path")
    end
    
    return fig
end

"""
Compute two-level system ground state energy using the actual SW results.

These are the coefficients from schrieffer_wolff(H, P; order=n) for
H = Δ/2 σz + g(σ+ + σ-), P = Subspace(σz => -1):

Order 2: H_P = -Δ/2 - g²/Δ
Order 4: H_P = -Δ/2 - g²/Δ + g⁴/Δ³
Order 6: H_P = -Δ/2 - g²/Δ + g⁴/Δ³ - 2g⁶/Δ⁵
Order 8: H_P = -Δ/2 - g²/Δ + g⁴/Δ³ - 2g⁶/Δ⁵ + 5g⁸/Δ⁷
Order 10: H_P = -Δ/2 - g²/Δ + g⁴/Δ³ - 2g⁶/Δ⁵ + 5g⁸/Δ⁷ - 14g¹⁰/Δ⁹
"""
function compute_two_level_energy(Δ::Float64, g::Float64, order::Int)
    E = -Δ/2
    
    if order >= 2
        E -= g^2 / Δ
    end
    if order >= 4
        E += g^4 / Δ^3
    end
    if order >= 6
        E -= 2*g^6 / Δ^5
    end
    if order >= 8
        E += 5*g^8 / Δ^7
    end
    if order >= 10
        E -= 14*g^10 / Δ^9
    end
    
    return E
end

"""
    figure4_rabi_bloch_siegert(; save_path=nothing)

Create figure comparing Jaynes-Cummings (RWA) with full Rabi model,
showing the Bloch-Siegert shift from counter-rotating terms.
"""
function figure4_rabi_bloch_siegert(; save_path=nothing)
    ω = 5.0   # Oscillator frequency
    Δ = 1.0   # Detuning
    g_range = range(0, 0.5, length=100)
    
    # JC (RWA) dispersive shift: χ_JC = -g²/Δ
    χ_JC = [-g^2/Δ for g in g_range]
    
    # Counter-rotating contribution: χ_CR = -g²/(Δ + 2ω)
    χ_CR = [-g^2/(Δ + 2*ω) for g in g_range]
    
    # Full Rabi (approximate): χ_full ≈ χ_JC + χ_CR
    χ_full = χ_JC .+ χ_CR
    
    fig = Figure(size = (800, 450))
    
    # Main plot
    ax1 = Axis(fig[1, 1],
        xlabel = L"g / \Delta",
        ylabel = L"\chi / \Delta",
        title = L"\mathrm{Dispersive shift: JC vs full Rabi } (\omega/\Delta = 5)",
    )
    
    lines!(ax1, g_range, χ_JC, color = COLORS[1], linewidth = 2.5,
           label = L"\mathrm{JC (RWA): } \chi = -g^2/\Delta")
    lines!(ax1, g_range, χ_full, color = COLORS[2], linewidth = 2.5,
           linestyle = :dash, label = L"\mathrm{Full Rabi: } \chi \approx -g^2/\Delta - g^2/(\Delta + 2\omega)")
    
    # Shade the Bloch-Siegert correction region
    band!(ax1, g_range, χ_JC, χ_full, color = (COLORS[3], 0.3))
    
    axislegend(ax1, position = :rt)
    
    # Inset: Bloch-Siegert shift magnitude
    ax2 = Axis(fig[1, 1],
        width = Relative(0.4),
        height = Relative(0.35),
        halign = 0.95,
        valign = 0.05,
        xlabel = L"g / \Delta",
        ylabel = "BS shift (%)",
        xlabelsize = 11,
        ylabelsize = 11,
        xticklabelsize = 9,
        yticklabelsize = 9,
        backgroundcolor = (:white, 0.9),
        title = "Bloch-Siegert correction",
        titlesize = 10,
    )
    
    bs_percent = 100 .* abs.(χ_CR) ./ abs.(χ_JC .+ 1e-10)  # Avoid div by 0
    lines!(ax2, g_range[2:end], bs_percent[2:end], color = COLORS[3], linewidth = 2)
    
    # For our parameters: BS shift is Δ/(Δ+2ω) = 1/11 ≈ 9.1%
    hlines!(ax2, [100*Δ/(Δ + 2*ω)], color = :gray, linestyle = :dash)
    
    if save_path !== nothing
        save(save_path, fig, px_per_unit = 3)
        println("Saved: $save_path")
    end
    
    return fig
end

"""
    figure5_combined_summary(; save_path=nothing)

Create a comprehensive summary figure with multiple panels showing
the power of higher-order SW transformations.
"""
function figure5_combined_summary(; save_path=nothing)
    fig = Figure(size = (1000, 800))
    
    Δ = 1.0
    
    # Panel (a): Two-level energy levels
    ax1 = Axis(fig[1, 1],
        xlabel = L"g / \Delta",
        ylabel = L"E / \Delta",
        title = "(a) Two-level system energy levels",
    )
    
    g_range = range(0, 0.6, length=80)
    E_g_exact = [-sqrt(Δ^2/4 + g^2) for g in g_range]
    E_e_exact = [sqrt(Δ^2/4 + g^2) for g in g_range]
    E_g_SW2 = [-Δ/2 - g^2/Δ for g in g_range]
    E_e_SW2 = [Δ/2 + g^2/Δ for g in g_range]
    E_g_SW4 = [-Δ/2 - g^2/Δ + g^4/Δ^3 for g in g_range]
    E_e_SW4 = [Δ/2 + g^2/Δ - g^4/Δ^3 for g in g_range]
    E_g_SW6 = [-Δ/2 - g^2/Δ + g^4/Δ^3 - 2*g^6/Δ^5 for g in g_range]
    E_e_SW6 = [Δ/2 + g^2/Δ - g^4/Δ^3 + 2*g^6/Δ^5 for g in g_range]
    
    lines!(ax1, g_range, E_g_exact, color = :black, linewidth = 2.5, label = "Exact")
    lines!(ax1, g_range, E_e_exact, color = :black, linewidth = 2.5)
    lines!(ax1, g_range, E_g_SW2, color = COLORS[1], linewidth = 2, linestyle = :dash, label = "Order 2")
    lines!(ax1, g_range, E_e_SW2, color = COLORS[1], linewidth = 2, linestyle = :dash)
    lines!(ax1, g_range, E_g_SW4, color = COLORS[2], linewidth = 2, linestyle = :dashdot, label = "Order 4")
    lines!(ax1, g_range, E_e_SW4, color = COLORS[2], linewidth = 2, linestyle = :dashdot)
    lines!(ax1, g_range, E_g_SW6, color = COLORS[3], linewidth = 2, linestyle = :dot, label = "Order 6")
    lines!(ax1, g_range, E_e_SW6, color = COLORS[3], linewidth = 2, linestyle = :dot)
    
    axislegend(ax1, position = :lt)
    
    # Panel (b): Energy gap
    ax2 = Axis(fig[1, 2],
        xlabel = L"g / \Delta",
        ylabel = L"(E_+ - E_-) / \Delta",
        title = "(b) Energy gap",
    )
    
    gap_exact = E_e_exact .- E_g_exact
    gap_SW2 = E_e_SW2 .- E_g_SW2
    gap_SW4 = E_e_SW4 .- E_g_SW4
    gap_SW6 = E_e_SW6 .- E_g_SW6
    
    lines!(ax2, g_range, gap_exact, color = :black, linewidth = 2.5, label = "Exact")
    lines!(ax2, g_range, gap_SW2, color = COLORS[1], linewidth = 2, linestyle = :dash, label = "Order 2")
    lines!(ax2, g_range, gap_SW4, color = COLORS[2], linewidth = 2, linestyle = :dashdot, label = "Order 4")
    lines!(ax2, g_range, gap_SW6, color = COLORS[3], linewidth = 2, linestyle = :dot, label = "Order 6")
    
    axislegend(ax2, position = :lt)
    
    # Panel (c): Relative error vs coupling
    ax3 = Axis(fig[2, 1],
        xlabel = L"g / \Delta",
        ylabel = "Relative error (%)",
        title = "(c) Approximation error in ground state",
    )
    
    err_SW2 = 100 .* abs.(E_g_SW2 .- E_g_exact) ./ abs.(E_g_exact)
    err_SW4 = 100 .* abs.(E_g_SW4 .- E_g_exact) ./ abs.(E_g_exact)
    err_SW6 = 100 .* abs.(E_g_SW6 .- E_g_exact) ./ abs.(E_g_exact)
    
    lines!(ax3, g_range, err_SW2, color = COLORS[1], linewidth = 2.5, label = "Order 2")
    lines!(ax3, g_range, err_SW4, color = COLORS[2], linewidth = 2.5, label = "Order 4")
    lines!(ax3, g_range, err_SW6, color = COLORS[3], linewidth = 2.5, label = "Order 6")
    
    axislegend(ax3, position = :lt)
    
    # Panel (d): JC dispersive shift and Kerr
    ax4 = Axis(fig[2, 2],
        xlabel = L"g / \Delta",
        ylabel = L"\mathrm{Coefficient} / \Delta",
        title = "(d) Jaynes-Cummings: dispersive shift and Kerr",
    )
    
    g_range_jc = range(0, 0.4, length=80)
    χ_SW2 = [-g^2/Δ for g in g_range_jc]
    χ_SW4 = [-g^2/Δ + 4*g^4/(3*Δ^3) for g in g_range_jc]
    K_SW4 = [4*g^4/(3*Δ^3) for g in g_range_jc]
    
    lines!(ax4, g_range_jc, χ_SW2, color = COLORS[1], linewidth = 2.5, 
           linestyle = :dash, label = L"\chi \mathrm{ (order 2)}")
    lines!(ax4, g_range_jc, χ_SW4, color = COLORS[2], linewidth = 2.5, 
           label = L"\chi \mathrm{ (order 4)}")
    lines!(ax4, g_range_jc, K_SW4, color = COLORS[3], linewidth = 2.5, 
           linestyle = :dot, label = L"K \mathrm{ (Kerr, order 4)}")
    
    hlines!(ax4, [0], color = :gray, linewidth = 1)
    axislegend(ax4, position = :lb)
    
    if save_path !== nothing
        save(save_path, fig, px_per_unit = 3)
        println("Saved: $save_path")
    end
    
    return fig
end

# Main execution
function generate_all_figures(output_dir::String = ".")
    set_publication_theme!()
    
    mkpath(output_dir)
    
    println("Generating publication figures for UnitaryTransformations.jl")
    println("="^60)
    
    println("\n1. Two-level system comparison...")
    figure1_two_level_system(save_path = joinpath(output_dir, "two_level_system.png"))
    
    println("\n2. Jaynes-Cummings dispersive regime...")
    figure2_jaynes_cummings(save_path = joinpath(output_dir, "jaynes_cummings.png"))
    
    println("\n3. Convergence analysis...")
    figure3_convergence(save_path = joinpath(output_dir, "convergence.png"))
    
    println("\n4. Rabi model / Bloch-Siegert shift...")
    figure4_rabi_bloch_siegert(save_path = joinpath(output_dir, "bloch_siegert.png"))
    
    println("\n5. Combined summary figure...")
    figure5_combined_summary(save_path = joinpath(output_dir, "summary.png"))
    
    # Also save as PDF for publication
    println("\n6. Saving PDF versions...")
    figure1_two_level_system(save_path = joinpath(output_dir, "two_level_system.pdf"))
    figure2_jaynes_cummings(save_path = joinpath(output_dir, "jaynes_cummings.pdf"))
    figure5_combined_summary(save_path = joinpath(output_dir, "summary.pdf"))
    
    println("\n" * "="^60)
    println("All figures generated successfully!")
    println("Output directory: $output_dir")
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    output_dir = length(ARGS) > 0 ? ARGS[1] : joinpath(@__DIR__)
    generate_all_figures(output_dir)
end
