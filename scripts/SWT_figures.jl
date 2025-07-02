using DrWatson
@quickactivate "2025-Williams-Skyrmions_frustrated_quantum_dimers"

using GLMakie
using Sunny, CairoMakie, JLD2
using CairoMakie
using Statistics

includet(srcdir("model_SWT.jl"))


# Load saved ground states
data = load("entangled_coherents_run_alpha06_60.jld2")
coherent_bufs = data["coherent_buffers"]

#############################################################
# Real space plot of individual example
#############################################################
GLMakie.activate!()
# Generate the parameters
N1 = 209
N2 = 28
Bs = range(0.0,1.0,N1)
αs = range(0.06,0.6,N2)
J=1.0
Δ=1.2
dims = (5,5,1)

αs[10]
Bs[161]
# Make a system
α = αs[25]
B = Bs[135]
    (; sys, crystal) = tl_dimer_model_from_low_energy_params(; dims, α, Δ) 
    set_field!(sys, [0, 0, B])
    coherents = coherent_bufs[135,25]
    for unit in Sunny.eachunit(sys)
        set_coherent!(sys, coherents[unit], unit)
    end
plot_spins(sys)


#############################################################
# Figure 5: YZ-S, SkX-I
#############################################################
idcs = [[181,3],[163,5]] # YZ-S, SkX-I

qs0 = [[0, 0, 0], [2/3, -1/3, 0], [1/2, 0, 0], [0, 0, 0]]
qs1 = [[0, 0, 1], [2/3, -1/3, 1], [1/2, 0, 1], [0, 0, 1]]
path0 = q_space_path(crystal, qs0, 400, labels = ["Γ", "K", "M", "Γ"])
path1 = q_space_path(crystal, qs1, 400, labels = ["Γ", "K", "M", "Γ"])
energies = range(0, 3.0, 400)

# Do SWT calculations
ress0 = []
ress1 = []
αs_plotted = Float64[]
Bs_plotted = Float64[]
for idx in idcs
    # Load up saved data
    α = αs[idx[2]]
    B = Bs[idx[1]]
    coherents = coherent_bufs[idx...]

    # Make a system
    (; sys, crystal) = tl_dimer_model_from_low_energy_params(; dims, α, Δ) 
    set_field!(sys, [0, 0, B])
    for unit in Sunny.eachunit(sys)
        set_coherent!(sys, coherents[unit], unit)
    end
    minimize_energy!(sys)

    # Make a spin wave theory
    measure = ssf_trace(sys)
    swt = SpinWaveTheory(sys; measure, regularization = 1e-6)

    # Calculate the intensities
    @time res0 = intensities(swt, path0; energies, kernel=gaussian(; fwhm=0.1))
    @time res1 = intensities(swt, path1; energies, kernel=gaussian(; fwhm=0.1))

    push!(ress0, res0)
    push!(ress1, res1)
    push!(αs_plotted, α)
    push!(Bs_plotted, B)
end


fig5 = with_theme(theme_latexfonts()) do
    fig = Figure(; size=(1200, 330))
    colorrange1 = (0, 150)
    colorrange2 = (0, 125)
    xticks = ress1[1].qpts.xticks
    qs = ress1[1].qpts.qs
    ylabel = "Energy (J)"
    xlabel = "Momentum"
    title0 = L"q_z = 0"
    title1 = L"q_z = 1"
    titlesize=24
    g1 = fig[1,1] = GridLayout()
    g2 = fig[1,2] = GridLayout()
    g3 = fig[1,3] = GridLayout()
    g4 = fig[1,4] = GridLayout()
    g5 = fig[1,5] = GridLayout()
    g6 = fig[1,6] = GridLayout()
    xticklabelsize=20
    yticklabelsize=20
    xlabelsize=22
    ylabelsize=22
    ax1 = Axis(g1[1,1]; xticks, ylabel, xlabel, title=title0, titlesize, yticklabelsize, xticklabelsize, xlabelsize, ylabelsize)
    ax2 = Axis(g2[1,1]; xticks, ylabel, xlabel, title=title1, titlesize, yticklabelsize, xticklabelsize, xlabelsize, ylabelsize)
    ax3 = Axis(g4[1,1]; xticks, ylabel, xlabel, title=title0, titlesize, yticklabelsize, xticklabelsize, xlabelsize, ylabelsize)
    ax4 = Axis(g5[1,1]; xticks, ylabel, xlabel, title=title1, titlesize, yticklabelsize, xticklabelsize, xlabelsize, ylabelsize)
    hideydecorations!(ax2)
    hideydecorations!(ax4)
    hm0 = heatmap!(ax1, 1:length(qs), energies, ress0[1].data'; colorrange=colorrange1, colormap=:gnuplot2)
    heatmap!(ax2, 1:length(qs), energies, ress1[1].data'; colorrange=colorrange1, colormap=:gnuplot2)
    hm1 = heatmap!(ax3, 1:length(qs), energies, ress0[2].data'; colorrange=colorrange2, colormap=:gnuplot2)
    heatmap!(ax4, 1:length(qs), energies, ress1[2].data'; colorrange=colorrange2, colormap=:gnuplot2)
    Colorbar(g3[1,1], hm0; ticklabelsize=20)
    Colorbar(g6[1,1], hm1; ticklabelsize=20)

    
    # for (label, layout) in zip(["a", "b", "c", "d"], [g1, g2, g4, g5])
    for (label, layout) in zip(["a", "b"], [g1, g4])

        Label(layout[1,1,TopLeft()], label,
        fontsize=22,
        font = :bold,
        padding = (0, 0, 10, 0),
        halign=:center)
    end

    fig   
end

wsave(plotsdir("PaperFigures", "Fig5.pdf"), fig5)

#########################################################################
# Figure 11: Single-E
#########################################################################
idcs = [[181,3],[163,5]] # YZ-S, SkX-I
grid = q_space_grid(crystal, [1, 0, 0], range(-0.4, 0.4, 200), [0, 1, 0], (-0.4, 0.4); offset = [0,0,1], orthogonalize=true)
E = 0.012

# Do SWT calculations
ress= []
αs_plotted = Float64[]
Bs_plotted = Float64[]
for idx in idcs
    # Load up saved data
    α = αs[idx[2]]
    B = Bs[idx[1]]
    coherents = coherent_bufs[idx...]

    # Make a system
    (; sys, crystal) = tl_dimer_model_from_low_energy_params(; dims, α, Δ) 
    set_field!(sys, [0, 0, B])
    for unit in Sunny.eachunit(sys)
        set_coherent!(sys, coherents[unit], unit)
    end
    minimize_energy!(sys)

    # Make a spin wave theory
    measure = ssf_trace(sys)
    swt = SpinWaveTheory(sys; measure, regularization = 1e-6)

    # Calculate the intensities
    @time res = intensities(swt, grid; energies=[E], kernel=gaussian(; fwhm=0.01))

    push!(ress, res)
    push!(αs_plotted, α)
    push!(Bs_plotted, B)
    
end
idx = idcs[1]
α = αs[idx[2]]
B = Bs[idx[1]]
coherents = coherent_bufs[idx...]
# Make a system
(; sys, crystal) = tl_dimer_model_from_low_energy_params(; dims, α, Δ) 
set_field!(sys, [0, 0, B])
for unit in Sunny.eachunit(sys)
    set_coherent!(sys, coherents[unit], unit)
end
minimize_energy!(sys)
GLMakie.activate!()
plot_spins(sys)
# Make a spin wave theory
measure = ssf_trace(sys)
swt = SpinWaveTheory(sys; measure, regularization = 1e-6)

rotations = [([0,0,1], n*(2π/3)) for n in 0:2]
weights = [1, 1, 1]
# Calculate the intensities
@time res = domain_average(crystal, grid; rotations, weights) do path_rotated
            intensities(swt, path_rotated; energies=[E], kernel=gaussian(; fwhm=0.01))
        end

push!(ress, res)
push!(αs_plotted, α)
push!(Bs_plotted, B)

CairoMakie.activate!()
fig11 = with_theme(theme_latexfonts()) do
    fig = Figure(size=(1500, 350)) #was (1000,350) for 2 panel
    g1 = fig[1,1] = GridLayout()
    g2 = fig[1,2] = GridLayout()
    #gblank = fig[1,2] = GridLayout()
    g3 = fig[1,3] = GridLayout()    
    axisopts = (;
        xticklabelsize=20,
        yticklabelsize=20,
        xlabelsize=22,
        ylabelsize=22,
        xlabel="[H, 0, 1]",
        ylabel="[-K/2, K, 1]",
    )
    plot_intensities!(g1, ress[1]; colorrange=(0, 8000), axisopts)
    plot_intensities!(g2, ress[4]; colorrange=(0, 3000), axisopts)
    plot_intensities!(g3, ress[2]; colorrange=(0, 3000), axisopts)
    for (label, layout) in zip(["a", "b", "c"], [g1, g2, g3])
        Label(layout[1,1,TopLeft()], label,
        fontsize=22,
        font = :bold,
        padding = (0, -120, 10, 0),
        halign=:center)
    end
    fig.content[2].ticklabelsize = 20
    fig.content[4].ticklabelsize = 20
    fig.content[6].ticklabelsize = 20
    #colsize!(fig.layout, 2, Auto(0.1))
    fig
end

wsave(plotsdir("PaperFigures", "Fig11v2.pdf"), fig11)

wsave(plotsdir("ress_for_Fig11.jld2"),"ress",ress)

#############################################################
# Figure 10: Remaining phases (Spiral phases and SkX-II)
#############################################################
# Fletcher-specified indices
#idcs = [[105,15] (AFM-CS), [163,5] (SkX-I), [181,3] (YZ-S), [129,13] (3Q-II), 
#[105,20] (2Q-II), [180,2] (2Q-I), [135,23] (FM-CS), [157,10] (SkX-II), [161,10] (SkX-II)]
#idcs = [[181,3],[105,15],[135,23],[163,5],[161,10]] #only single Q and SkX's
#idcs = [[43,5],[43,7],[43,9],[43,10], [43,11],[43,12],[43,13]] #transition from QPM to AFM-CS to see singlet-triplon dispersion
#idcs = [[1,5],[1,7],[1,9],[1,10], [1,11],[1,12],[1,13]] #transition from QPM to AFM-CS to see singlet-triplon dispersion
#idcs = [[1,15],[1,17],[1,19]]
#[181,3] = YZ-Spiral, [1,18] = zero field AFM-CS, 
# idcs = [[181,3],[163,5]]
#idcs =[[181,3]] = YZ-S, #[163,5] = in the SkX-I phase, [161,10] = in the SkX-II phase

idcs = [
    [105, 15], # AFM-CS
    [135,23],  # FM-CS
    [129,13],  # 3Q-II
    [157,10],  # SkX-II
]
αs[13]
Bs[129]
qs0 = [[0, 0, 0], [2/3, -1/3, 0], [1/2, 0, 0], [0, 0, 0]]
qs1 = [[0, 0, 1], [2/3, -1/3, 1], [1/2, 0, 1], [0, 0, 1]]
path0 = q_space_path(crystal, qs0, 300, labels = ["Γ", "K", "M", "Γ"])
path1 = q_space_path(crystal, qs1, 300, labels = ["Γ", "K", "M", "Γ"])
energies_all = [
    range(0, 3.0, 200),
    range(0, 7.0, 200),
    range(0, 4.0, 200),
    range(0, 4.0, 200),
]

# Do SWT calculations
ress0 = []
ress1 = []
αs_plotted = Float64[]
Bs_plotted = Float64[]
for (n, idx) in enumerate(idcs)
    # Load up saved data
    α = αs[idx[2]]
    B = Bs[idx[1]]
    coherents = coherent_bufs[idx...]

    # Make a system
    (; sys, crystal) = tl_dimer_model_from_low_energy_params(; dims, α, Δ) 
    set_field!(sys, [0, 0, B])
    for unit in Sunny.eachunit(sys)
        set_coherent!(sys, coherents[unit], unit)
    end
    minimize_energy!(sys)

    # Make a spin wave theory
    measure = ssf_trace(sys)
    swt = SpinWaveTheory(sys; measure, regularization = 1e-6)

    # Calculate the intensities
    @time res0 = intensities(swt, path0; energies=energies_all[n], kernel=gaussian(; fwhm=0.1))
    @time res1 = intensities(swt, path1; energies=energies_all[n], kernel=gaussian(; fwhm=0.1))

    push!(ress0, res0)
    push!(ress1, res1)
    push!(αs_plotted, α)
    push!(Bs_plotted, B)
end


fig10 = with_theme(theme_latexfonts()) do
    fig = Figure(; size=(1200, 660))

    g1 = fig[1,1] = GridLayout()
    g2 = fig[1,2] = GridLayout()
    g3 = fig[1,3] = GridLayout()

    g4 = fig[1,4] = GridLayout()
    g5 = fig[1,5] = GridLayout()
    g6 = fig[1,6] = GridLayout()

    g7 = fig[2,1] = GridLayout()
    g8 = fig[2,2] = GridLayout()
    g9 = fig[2,3] = GridLayout()

    g10 = fig[2,4] = GridLayout()
    g11 = fig[2,5] = GridLayout()
    g12 = fig[2,6] = GridLayout()

    colorrange1 = (0, 150)
    colorrange2 = (0, 150)
    colorrange3 = (0, 125)
    colorrange4 = (0, 125)
    xticks = ress1[1].qpts.xticks
    qs = ress1[1].qpts.qs
    ylabel = "Energy (J)"
    xlabel = "Momentum"
    title0 = L"q_z = 0"
    title1 = L"q_z = 1"
    titlesize=24
    xticklabelsize=20
    yticklabelsize=20
    xlabelsize=22
    ylabelsize=22

    ax1 = Axis(g1[1,1]; yticks=0:3, xticks, ylabel, xlabel, title=title0, titlesize, yticklabelsize, xticklabelsize, xlabelsize, ylabelsize)
    ax2 = Axis(g2[1,1]; xticks, ylabel, xlabel, title=title1, titlesize, yticklabelsize, xticklabelsize, xlabelsize, ylabelsize)
    ax3 = Axis(g4[1,1]; yticks=0:7, xticks, ylabel, xlabel, title=title0, titlesize, yticklabelsize, xticklabelsize, xlabelsize, ylabelsize)
    ax4 = Axis(g5[1,1]; xticks, ylabel, xlabel, title=title1, titlesize, yticklabelsize, xticklabelsize, xlabelsize, ylabelsize)
    ax5 = Axis(g7[1,1]; yticks=0:4, xticks, ylabel, xlabel, title=title0, titlesize, yticklabelsize, xticklabelsize, xlabelsize, ylabelsize)
    ax6 = Axis(g8[1,1]; xticks, ylabel, xlabel, title=title1, titlesize, yticklabelsize, xticklabelsize, xlabelsize, ylabelsize)
    ax7 = Axis(g10[1,1]; yticks=0:4, xticks, ylabel, xlabel, title=title0, titlesize, yticklabelsize, xticklabelsize, xlabelsize, ylabelsize)
    ax8 = Axis(g11[1,1]; xticks, ylabel, xlabel, title=title1, titlesize, yticklabelsize, xticklabelsize, xlabelsize, ylabelsize)

    hideydecorations!(ax2)
    hideydecorations!(ax4)
    hideydecorations!(ax6)
    hideydecorations!(ax8)

    hm0 = heatmap!(ax1, 1:length(qs), energies_all[1], ress0[1].data'; colorrange=colorrange1, colormap=:gnuplot2)
    heatmap!(ax2, 1:length(qs), energies_all[1], ress1[1].data'; colorrange=colorrange1, colormap=:gnuplot2)

    hm1 = heatmap!(ax3, 1:length(qs), energies_all[2], ress0[2].data'; colorrange=colorrange2, colormap=:gnuplot2)
    heatmap!(ax4, 1:length(qs), energies_all[2], ress1[2].data'; colorrange=colorrange2, colormap=:gnuplot2)

    hm2 = heatmap!(ax5, 1:length(qs), energies_all[3], ress0[3].data'; colorrange=colorrange3, colormap=:gnuplot2)
    heatmap!(ax6, 1:length(qs), energies_all[3], ress1[3].data'; colorrange=colorrange3, colormap=:gnuplot2)

    hm3 = heatmap!(ax7, 1:length(qs), energies_all[4], ress0[4].data'; colorrange=colorrange4, colormap=:gnuplot2)
    heatmap!(ax8, 1:length(qs), energies_all[4], ress1[4].data'; colorrange=colorrange4, colormap=:gnuplot2)

    Colorbar(g3[1,1], hm0; ticklabelsize=20)
    Colorbar(g6[1,1], hm1; ticklabelsize=20)
    Colorbar(g9[1,1], hm2; ticklabelsize=20)
    Colorbar(g12[1,1], hm3; ticklabelsize=20)

    
    for (label, layout) in zip(["a", "b", "c", "d"], [g1, g4, g7, g10])
        Label(layout[1,1,TopLeft()], label,
        fontsize=22,
        font = :bold,
        padding = (0, 0, 10, 0),
        halign=:center)
    end
    fig   
end

wsave(plotsdir("PaperFigures", "Fig10v2.pdf"), fig10)