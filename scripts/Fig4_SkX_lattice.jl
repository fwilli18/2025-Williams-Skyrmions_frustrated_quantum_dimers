using DrWatson
@quickactivate "2025-Williams-Skyrmions_frustrated_quantum_dimers"
using GLMakie, BSplineKit, Sunny, LinearAlgebra, JLD2, Revise, CairoMakie
GLMakie.activate!()

includet(srcdir("model.jl"))
includet(srcdir("helpers.jl"))
includet(srcdir("plotting.jl"))
# filename =  "C:\\Users\\Fletcher\\Documents\\Batista\\Coupled Dimers\\TriangularLatticeDimers\\entangled_coherents_run_alpha06_60.jld2"

#saved_coherents = load(datadir("entangled_coherents_run_alpha06_60.jld2"),"coherent_buffers")
saved_coherents = load("entangled_coherents_run_alpha06_60.jld2","coherent_buffers")


N1 = 209
N2 = 28
Bs = range(0.0,1.0,N1)
αs = range(0.06,0.6,N2)
J=1.0
Δ=1.2
dims = (5,5,1)

GLMakie.activate!()
CairoMakie.activate!()
# Pick one of the systems
B_index = 161
alpha_index = 10
Bs[B_index]
αs[alpha_index]

sys = entangled_model_from_alpha_b(; dims, J, Δ, α=αs[alpha_index], B = Bs[B_index])
for site in Sunny.eachsite(sys[1].sys)
    set_coherent!(sys[1], saved_coherents[B_index, alpha_index][site], site)
end
sys_big = repeat_periodically(sys[1], (5, 5, 1))
#sys_big2 = repeat_periodically(sys, (3, 3, 1))


#Spin configuration with chirality figure
fig = let
    fig = Figure(; size=(900, 600))
    ax = LScene(fig[1,1], show_axis=false) 
    plot_triangular_plaquettes!(ax, sun_berry_curvature, sys_big, 9, 18, 8, 17; colorrange=(-0.5, 0.5), colormap=(:PRGn, 0.5))
    sys_layer = layer_from_bilayer(sys_big,2)
    plot_spins_alt!(ax, sys_layer; color = [s[3] for s in sys_layer.dipoles],
                                #colorfn = i -> norm(sys_layer.dipoles[i]), 
                                colorrange=(0, 0.5), show_cell=false, ndims=2, camdist=5.0)
    Colorbar(fig[1, 2], limits = (0, 0.5), colormap = :viridis,
    #Colorbar(fig[2, 1], limits = (0, 0.5), colormap = :viridis,
        flipaxis = true,
        #flipaxis = false, vertical = false,
        label=L"S^z", labelsize =40, size = 10, ticklabelsize =25)
    Colorbar(fig[2, 1], limits = (-1.5, 1.5), colormap = (:PRGn, 0.5),
        vertical = false, label=L"ρ_{jkl}", labelsize =40, ticklabelsize =25, flipaxis = false, size = 10)
    
    fig
end

#save("SkXII-top-v2.png",fig)
#SWT
measure = ssf_perp(sys[1])
swt = SpinWaveTheory(sys[1]; measure)
qs = [[-1/2, 0, 0], [0, 0, 0], [1/2, 1/2, 0]]
path = q_space_path(sys[2], qs, 400)
res = intensities_bands(swt, path)
plot_intensities(res)



itemp0, atemp0, btemp0 = structure_factor_2D(sys_big.sys; c =0.0,
                alims=(-2/3, 2/3), blims=(-1/√3, 1/√3), npoints=500,contraction=:full)
        
itemp5, atemp5, btemp5 = structure_factor_2D(sys_big; c =0.51,
                alims=(-2/3, 2/3), blims=(-1/√3, 1/√3), npoints=500,contraction=:full)


            
GLMakie.activate!()
CairoMakie.activate!()
#Single SSF plot for either z component or transverse components
begin
    key_index = 2
    
    cmin = minimum([minimum(keysfcomps[key_index,1]), minimum(keysfcomps[key_index,2]), minimum(keysfcomps[key_index,3]),
                       minimum(keysfcomps[key_index+11,1]), minimum(keysfcomps[key_index+11,2]), minimum(keysfcomps[key_index+11,3])])
    cmax = 4
    #this for S_perp
    fig_p = Figure(backgroundcolor=:transparent)
    qzp=1 #Pick 0 or 1
    title_fig = L"S^{\perp}(q_x , q_y)"

    ax_p = Axis(fig_p[1,1]; title=title_fig,
        #xlabel = L"q_{x}",
        #ylabel = L"q_{y}",
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        xlabelsize = 36,
        ylabelsize = 36,
        xticklabelsize = 30,
        yticklabelsize = 30,
        limits = (-0.5,0.5,-0.5,0.5),
        titlesize = 42,#"qz = 0, α=$(round(αs[alpha_index],digits=4)), B=$(round(Bs[B_index], digits=4))", 
        aspect=true)
    hm_p = heatmap!(ax_p, keyas[key_index], keybs[key_index], 0.5*(keysfcomps[key_index+11*qzp,1]+keysfcomps[key_index+11*qzp,2]),
                    colormap=:Blues,colorrange=(cmin,cmax))
    #Colorbar(fig_p[:,end+1],hm_p)
    colsize!(fig_p.layout, 1, Aspect(1, 0.6183))
    text!(ax_p, -0.22, -0.47; text=L"q_z = 1", fontsize=40)
    
    
    #this for S^{zz}
    fig_z = Figure(backgroundcolor = :transparent)
    qzz=0 #Pick 0 or 1
    title_fig = L"S^{zz}(q_x , q_y)"
    ax_z = Axis(fig_z[1,1]; title=title_fig,
        #xlabel = L"q_{x}",
        #ylabel = L"q_{y}",
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        xlabelsize = 36,
        ylabelsize = 36,
        xticklabelsize = 30,
        yticklabelsize = 30,
        titlesize = 42, 
        aspect=true)
    hm_z = heatmap!(ax_z, keyas[key_index], keybs[key_index], keysfcomps[key_index+11*qzz,3],colormap=:Blues, colorrange=(cmin,cmax))
    #Colorbar(fig_z[:,end+1],hm_z)
    colsize!(fig_z.layout, 1, Aspect(1, 0.6183))
    #hidespines!(ax_p)
    text!(ax_z, -0.28, -0.47; text=L"q_z = 0", fontsize=40)

end

    fig_p
    fig_z

    save("Skx2_ssf_perp.png",fig_p)
    save("Skx2_ssf_zz.png",fig_z)