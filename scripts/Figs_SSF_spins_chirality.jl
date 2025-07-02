using DrWatson
@quickactivate "2025-Williams-Skyrmions_frustrated_quantum_dimers"
using GLMakie, BSplineKit, Sunny, LinearAlgebra, JLD2, Revise
#using CairoMakie

includet(srcdir("model.jl"))
includet(srcdir("helpers.jl"))

#filename =  "C:\\Users\\Fletcher\\Documents\\Batista\\Coupled Dimers\\TriangularLatticeDimers\\data\\entangled_coherents_run_alpha06_60.jld2"
#saved_coherents = load(filename,"coherent_buffers")
saved_coherents = load("entangled_coherents_run_alpha06_60.jld2","coherent_buffers")


begin
    N1 = 209
    N2 = 28
    Bs = range(0.0,1.0,N1)
    αs = range(0.06,0.6,N2)
    J=1.0
    Δ=1.2
    dims = (5,5,1)
    dims_full = (5,5,1,2)
    dims_ent = (5,5,1,1)
    energies = zeros(N1, N2)
end


#Allocate and initalize the systems and corresponding energies
#begin
    syss_ent = [entangled_model_from_alpha_b(; dims, J, Δ, α=αs[j], B = Bs[i]) for i in 1:N1, j in 1:N2]
    for j in 1:N2, i in 1:N1, site in Sunny.eachsite(syss_ent[i,j][1].sys)
        set_coherent!(syss_ent[i,j][1],saved_coherents[i,j][site],site)
    end
    bigger_syss_ent = [repeat_periodically(syss_ent[i,j][1],(2,2,1)) for i in 1:N1, j in 1:N2]
#end


#Compute singlet character over various sites within skyrmion region
chars = []
temp = 0
P = 0.5*[2 0 0 0; 0 1 -1 0; 0 -1 1 0; 0 0 0 0]
#P is the projection matrix onto the {|1,1>, |0,0>} subspace
Sk_points = [(187,3),(173,6),(161,9),(147,13),(139,16),(133,18),(127,20)]
length(Sk_points)
begin
    for i = 1:length(Sk_points)
        B_ind = Sk_points[i][1]
        alph_ind = Sk_points[i][2]
        for site in Sunny.eachsite(syss_ent[B_ind,alph_ind].sys)
            psi = saved_coherents[B_ind,alph_ind][site]
            temp = temp + psi' * P * psi
        end
        push!(chars, temp/25)
        temp = 0
    end
end


#Plot spins

    B_index = 127
    alpha_index = 20
    fig = Figure()
    ax = LScene(fig[1,1], show_axis=false) #; title="Δ=$(round(Δs[n], digits=4))")
    Sunny.Plotting.plot_spins!(ax, bigger_syss_ent[B_index,alpha_index].sys_origin; 
                                colorfn = i -> norm(bigger_syss_ent[B_index,alpha_index].sys_origin.dipoles[i]), 
                                colorrange=(0, 0.5))
    Colorbar(fig[1, 2], limits = (0, 0.5), colormap = :viridis)
    fig


psi = saved_coherents[B_index,alpha_index][2,2]
psi' * P * psi


Ms = Array{Float64}(undef,N1,N2) 
npoints = 500
is = Array{Matrix{Float64}}(undef,N1,N2)
as = Array{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}(undef,N1,N2)
bs = Array{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}(undef,N1,N2)
#Structure factors and magnetizations
for i in 1:N1, j in 1:N2
    Ms[i,j] = norm(magnetization(syss_ent[i,j].sys_origin))
    #is[i,j], as[i,j], bs[i,j] = structure_factor_2D(bigger_syss_ent[i,j]; c =0.51,
    #alims=(-2/3, 2/3), blims=(-1/√3, 1/√3), npoints=500)
end


#Locations of key points in phase Diagram

#Before change:keypoints = [(105,15),(163,5),(181,3),(129,13),(105,20),(180,2),(135,23),(157,10),(161,10),(121,13),(111,20)]
#After change:
keypoints = [(105,15),(108,20),(181,3),(129,13),(105,20),(181,2),(135,23),(101,17),(105,17),(121,13),(111,20)]
    keyis = []
    keyas = []
    keybs = [] 
    for z = 1:2
        for point in keypoints
            if z==1
                itemp, atemp, btemp = structure_factor_2D(bigger_syss_ent[point[1],point[2]]; c =0.0,
                alims=(-2/3, 2/3), blims=(-1/√3, 1/√3), npoints=500,contraction=:full)
            else
                itemp, atemp, btemp = structure_factor_2D(bigger_syss_ent[point[1],point[2]]; c =0.51,
                alims=(-2/3, 2/3), blims=(-1/√3, 1/√3), npoints=500,contraction=:full)
            end
            push!(keyis, itemp)
            push!(keyas, atemp)
            push!(keybs, btemp)
        end
    end

keyixs = []
keyiys = []
keyizs = []

for (n, intens) in enumerate(keyis)
    xtemp = Array{Float64}(undef,500,500)
    ytemp = Array{Float64}(undef,500,500)
    ztemp = Array{Float64}(undef,500,500)
    for j in 1:500, k in 1:500
        xtemp[j,k] = intens[j,k][1,1]
        ytemp[j,k] = intens[j,k][2,2]
        ztemp[j,k] = intens[j,k][3,3]
    end
    push!(keyixs,xtemp)
    push!(keyiys,ytemp)
    push!(keyizs,ztemp)

end

keysfcomps = Array{Matrix{Float64}}(undef,length(keyixs),3)
for i = 1:length(keyixs)
    keysfcomps[i,1] = keyixs[i]
    keysfcomps[i,2] = keyiys[i]
    keysfcomps[i,3] = keyizs[i]
end
keysyss = []
for i = 1:length(keyixs)
    push!(keysyss, bigger_syss_ent)
end

Ms = load("Ms_209_28_opt_June18alphas06_30.jld2","Ms")
is = load("sfs_209_28_opt_June18alphas06_30.jld2","is")
as = load("as_209_28_opt_June18alphas06_30.jld2","as")
bs = load("bs_209_28_opt_June18alphas06_30.jld2","bs")



#Single magnetization curve
begin
    f = Figure()
    numcols = 1
    α_start = 2
    for α_ind in 1:1
        ax = Axis(f[fldmod1(α_ind, numcols)...]; title="M vs B (α=$(αs[α_ind-1+α_start]))", aspect=true)
        lines!(ax, Bs[:], Ms[:,(α_ind-1+α_start)], color =:blue)
    end
    f
end

#Magnetization Curves
begin
    f = Figure()
    numcols = 5
    α_start = 10
    for α_ind in 1:10
        ax = Axis(f[fldmod1(α_ind, numcols)...]; title="M vs B (α=$(αs[α_ind-1+α_start]))", aspect=true)
        lines!(ax, Bs[:], Ms[:,(α_ind-1+α_start)], color =:blue)
    end
    f
end
    save("M_B_curve_10alphas.png",f)

#Plot spins

    B_index = 181
    alpha_index = 3
    fig = Figure()
    ax = LScene(fig[1,1], show_axis=false) #; title="Δ=$(round(Δs[n], digits=4))")
    Sunny.plot_spins!(ax, bigger_syss_ent[B_index,alpha_index].sys_origin;
                                color = [s[3] for s in bigger_syss_ent[B_index,alpha_index].sys_origin.dipoles],
                                #colorfn = i -> norm(bigger_syss_ent[B_index,alpha_index].sys_origin.dipoles[i]), 
                                colorrange=(0, 0.5))
    Colorbar(fig[1, 2], limits = (0, 0.5), colormap = :viridis)
    fig

# Plot structure factors
#begin
    fig = Figure()
    numcols = 5
    B_index = 170
    alpha_index = 4
    D_index = 50
    for (n, sf) in enumerate(is[B_index:B_index+14,alpha_index])
    #for (n, sf) in enumerate(is[B_index:B_index+14,D_index])
        ax = Axis(fig[fldmod1(n, numcols)...]; title="α=$(round(αs[alpha_index],digits=4)), B=$(round(Bs[B_index-1+n], digits=4))", aspect=true)
        #ax = Axis(fig[fldmod1(n, numcols)...]; title="Δ=$(round(Δs[D_index],digits=4)), B=$(round(Bs[B_index-1+n], digits=4))", aspect=true)
        
        heatmap!(ax, as[1,1], bs[1,1], sf)
    end
    fig
#end


# Plot berry curvature
text = "B = $(B_index), α=$(D_index)"
scene = plot_berry_curvature(bigger_syss_ent[B_index,alpha_index]; text, offset_spacing = 2)
αs[9]
begin
    B_ind = 161
    alpha_index = 10
    #D_index = 70
    texts = ["B = $(Bs[i]), α=$(αs[alpha_index])" for i in B_ind:(B_ind+14)]
    scene = plot_berry_curvature(bigger_syss_ent[B_ind:(B_ind+14),alpha_index]; size=(1200, 1200), numcols=5, texts, offset_spacing = 2)
    #texts = ["B = $(Bs[i]), Δ=$(Δs[D_index])" for i in B_ind:(B_ind+14)]
    #scene = plot_berry_curvature(bigger_syss_ent[B_ind:(B_ind+14),D_index]; size=(1200, 1200), numcols=5, texts, offset_spacing = 2)

end
#save("BerryCurve_alpha0_20Bs101.png",scene)


#Single SSF plot with color bar
begin
    fig = Figure()
    B_index = 161
    alpha_index = 10
    ax = Axis(fig[1,1]; title="qz = π, α=$(round(αs[alpha_index],digits=4)), B=$(round(Bs[B_index], digits=4))", aspect=true)
    hm = heatmap!(ax, as[1,1], bs[1,1], is[B_index,alpha_index])
    Colorbar(fig[:,end+1],hm)
    colsize!(fig.layout, 1, Aspect(1, 0.6183))
    fig
end
#save("sfs_alpha0_3B61.png",fig)

#Single Berry curvature

#begin
    B_ind = 181
    alpha_index = 2
    text = ["B = $(B_index), α=$(alpha_index)"]
    scene = plot_berry_curvature(bigger_syss_ent[B_ind,alpha_index]; size = (1200,1200), numcols=1, text, offset_spacing = 4)
    #scene = plot_berry_curvature(bigger_syss_ent[B_ind:(B_ind+14),alpha_index]; size=(1200, 1200), numcols=5, texts, offset_spacing = 2)
#end
#save("BerryCurve_alpha44B51.png",scene)


using CairoMakie
using GLMakie
GLMakie.activate!()
CairoMakie.activate!()
#Single SSF plot for either z component or transverse components
begin
    key_index = 1
    
    cmin = minimum([minimum(keysfcomps[key_index,1]), minimum(keysfcomps[key_index,2]), minimum(keysfcomps[key_index,3]),
                       minimum(keysfcomps[key_index+11,1]), minimum(keysfcomps[key_index+11,2]), minimum(keysfcomps[key_index+11,3])])
    #cmax = maximum([maximum(keysfcomps[key_index,1]), maximum(keysfcomps[key_index,2]), maximum(keysfcomps[key_index,3]),
    #                   maximum(keysfcomps[key_index+11,1]), maximum(keysfcomps[key_index+11,2]), maximum(keysfcomps[key_index+11,3])])
    cmax = 10
    #this for S_perp
    fig_p = Figure(backgroundcolor=:transparent)#Figure()#
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
        backgroundcolor=:white, 
        aspect=true)
    hm_p = heatmap!(ax_p, keyas[key_index], keybs[key_index], 0.5*(keysfcomps[key_index+11*qzp,1]+keysfcomps[key_index+11*qzp,2]),
                    colormap=:Blues,colorrange=(cmin,cmax))
    #Colorbar(fig_p[:,end+1],hm_p)
    colsize!(fig_p.layout, 1, Aspect(1, 0.6183))
    text!(ax_p, -0.22, -0.47; text=L"q_z = 1", fontsize=40)
    
    
    #this for S^{zz}
    fig_z = Figure(backgroundcolor = :transparent)#Figure()##
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
        backgroundcolor=:white,
        aspect=true)
    hm_z = heatmap!(ax_z, keyas[key_index], keybs[key_index], keysfcomps[key_index+11*qzz,3],colormap=:Blues, colorrange=(cmin,cmax))
    #Colorbar(fig_z[:,end+1],hm_z)
    colsize!(fig_z.layout, 1, Aspect(1, 0.6183))
    #hidespines!(ax_p)
    text!(ax_z, -0.28, -0.47; text=L"q_z = 0", fontsize=40)

end

    fig_p
    fig_z

    

    save("AFM-CS_ssf_perp.png",fig_p)
    save("AFM-CS_ssf_zz.png",fig_z)
#6 SSF plot with color bar for each coordinate
titles = [L"S^{xx}(q_x,q_y)", L"S^{yy}(q_x,q_y)", L"S^{zz}(q_x,q_y)"]

    cmin = minimum([minimum(keysfcomps[key_index,1]), minimum(keysfcomps[key_index,2]), minimum(keysfcomps[key_index,3]),
                       minimum(keysfcomps[key_index+11,1]), minimum(keysfcomps[key_index+11,2]), minimum(keysfcomps[key_index+11,3])])
    cmax = 5#maximum([maximum(keysfcomps[key_index,1]), maximum(keysfcomps[key_index,2]), maximum(keysfcomps[key_index,3]),
            #           maximum(keysfcomps[key_index+9,1]), maximum(keysfcomps[key_index+9,2]), maximum(keysfcomps[key_index+9,3])])
    for i = 1:3
        ax = Axis(fig[1,i];#fig[1,i]; 
            title=titles[i],#L"S^{zz}",#
            titlesize=20,
            aspect=true)
            if i!=3
                heatmap!(ax, as[1,1], bs[1,1], keysfcomps[key_index,i], colorrange=(cmin,cmax))#,colormap=:balance)
            else
                hm = heatmap!(ax, as[1,1], bs[1,1], keysfcomps[key_index,i], colorrange=(cmin,cmax))
            end

    end
    
    for i = 1:3
        ax2 = Axis(fig[2,i]; 
            title=titles[i],#L"S^{xx}+S^{yy}",
            titlesize=20,
            aspect=true)
        hm = heatmap!(ax2, as[1,1], bs[1,1], keysfcomps[key_index+11,i], colorrange=(cmin,cmax))
    end
    
    Colorbar(fig[:,end+1],hm, colorrange = (cmin,cmax))
    supertitle = Label(fig[0, :],
        L"α=%$(round(αs[keypoints[key_index][2]],digits=4)), B=%$(round(Bs[keypoints[key_index][1]], digits=4))",
        fontsize = 24)
    sideinfo = Label(fig[1, 0], L"q_z=0", rotation = pi/2, fontsize=20)
    sideinfo = Label(fig[2, 0], L"q_z=\pi", rotation = pi/2,  fontsize=20)
    #colsize!(fig.layout, 1, Aspect(1, 0.6183))
    rowsize!(fig.layout, 1, Aspect(2,1))
    rowsize!(fig.layout, 2, Aspect(2,1))
    fig



save("sfcomps_alpha44B51.png",fig)

begin
    plot_triangular_plaquettes!(ax, sun_berry_curvature, sys_big, 6,15 , 5, 14; colorrange=(-0.5,0.5))
    Colorbar(f[1, 2], limits = (0.0, 0.5), #colormap = :default,
        flipaxis = true, label=L"S^z", labelsize =40, size = 15, ticklabelsize =25,)
    Colorbar(f[2, 1], limits = (-1.5, 1.5), colormap = (:PRGn, 0.5),
        vertical = false, label=L"\chi_{ijk}", labelsize =40, ticklabelsize =25, flipaxis = false, size = 15)
    hidespines!(ax)
    hidedecorations!(ax)
    f
end
save("skx1_spins_chirality_bottom_v3.png",f)

