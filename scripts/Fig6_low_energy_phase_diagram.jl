using DrWatson
using CairoMakie
using GLMakie, BSplineKit, Sunny, LinearAlgebra, JLD2
@quickactivate "2025-Williams-Skyrmions_frustrated_quantum_dimers"

includet(srcdir("model.jl"))
includet(srcdir("helpers.jl"))
filename =  "dipoles_from_2200_pt_PS.jld2"
saved_dipoles = load(filename, "dipole_buffers")

begin
    N1 = 55
    N2 = 40
    hs = range(0,0.55,N1)
    Δs = range(0.95,1.35,N2)
    dims = (5,5,1)
    energies = zeros(N1, N2)
end

#Allocate and initalize the systems and corresponding energies
begin
    J1 = -1.0
    J2 = 2/(1+sqrt(5))
    #Initialize syss and allocate at same time
    syss = [triangular_effective_model_nnn(; dims, J1, J2, Δ = Δs[j], h = hs[i])[1] for i in 1:N1, j in 1:N2]
    for j in 1:N2, i in 1:N1, site in Sunny.eachsite(syss[i,j])
            set_dipole!(syss[i,j],saved_dipoles[i,j][site],site)
    end
    
end

#Calculate the magnetization of each system
Ms = Array{Float64}(undef,N1,N2) 
for i in 1:N1, j in 1:N2
    Ms[i,j] = norm(magnetization(syss[i,j]))
end

#If you want to view the static structure factors, uncomment out the next lines of code
#=npoints = 500
sfs = Array{Matrix{Float64}}(undef,N1,N2) 
as = Array{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}(undef,N1,N2)
bs = Array{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}(undef,N1,N2)
c = 1.0
for i in 1:N1, j in 1:N2
    sfs[i,j], as[i,j], bs[i,j] = structure_factor_2D(bigger_syss[i,j]; alims=(-2c/3, 2c/3), blims=(-1c/√3, 1c/√3), npoints=500)
end
# Plot structure factors
begin
    fig = Figure()
    numcols = 5
    h_index = 15
    del_index = 25
    cmin = minimum(sfs[h_index,del_index])
    cmax = maximum(sfs[h_index,del_index])
    
    for (n, sf) in enumerate(sfs[h_index:h_index+10,del_index])
        ax = Axis(fig[fldmod1(n, numcols)...]; title="K=$(round(Δs[del_index],digits=4)), h=$(round(hs[h_index-1+n], digits=4))", aspect=true)
        if minimum(sf) < cmin
            cmin = minimum(sf) 
        end
        if maximum(sf) > cmax
            cmax = maximum(sf)
        end
        colorrange = (cmin, 0.25*cmax)
        heatmap!(ax, sf; colorrange)
        
    end
    
    #Colorbar(fig[:,end+1],hm, ticklabelsize = 24)
    fig
end
=#



boundary_pts = [(0.95,0.0)] #Will store candidate positions of phase boundaries
slope_buffer = (Ms[2,1]-Ms[1,1])/(hs[2]-hs[1]) #First slope of magnetization-field curve
#determine locations of 1st order phase transitions
for j in 1:N2
    for i in 1:(N1-1)
        slope_next = (Ms[i+1,j]-Ms[i,j])/(hs[i+1]-hs[i])
        #println("The slope between $i and $(i+1) is $slope_next.")
        #check if slope of next pair changes by more than 30% of slope of pair before
        #check to make sure not fully polarized
        if abs((slope_next-slope_buffer)/slope_buffer) > 0.3 && (Ms[i,j]<0.99999999999 && Ms[i+1,j]<1.0)
            push!(boundary_pts, (Δs[j],hs[i]))
            #println("slope_buffer = $slope_buffer")
            #println("slope_next = $slope_next")
        end
        slope_buffer = slope_next
    end
    slope_buffer = (Ms[2,j]-Ms[1,j])/(hs[2]-hs[1])
end

#Add 2nd order transition points by inspection of structure factors
sec_order_pts = [(1.0,0.43),(Δs[7],hs[47]), (Δs[8],hs[48]), (Δs[9],hs[48]), (Δs[11],hs[46]),
(Δs[12],hs[45]), (Δs[13],hs[44]), (Δs[14],hs[43]), (Δs[15],hs[42]), (Δs[16],hs[41]),#boundary_pts[76], (Δs[17],hs[40]),#(Δs[10],hs[46]), 
(Δs[18],hs[39])]
for i in eachindex(sec_order_pts)
    push!(boundary_pts, sec_order_pts[i])
end


# Visualize candidate locations for phase boundaries
begin
    fig = Figure()
    ax = Axis(fig[1,1];
    xlabel = "Δ",
    ylabel = "h",
    xlabelsize = 20,
    ylabelsize = 20,
    xticklabelsize = 16,
    yticklabelsize = 16,
    title = "Coupled Dimers Low Energy Phase Diagram",
    titlesize = 22,
    limits = (0.95,1.35,0.0,0.6))
    scatter!(ax, boundary_pts)
    fig
end

#Collect boundary points for interpolations
#First make the VS/M boundary
vs_m_pts_indices = [6,39,47,106]
vs_m_pts = [boundary_pts[2]]
for i in eachindex(vs_m_pts_indices)
    push!(vs_m_pts, boundary_pts[vs_m_pts_indices[i]])
end
vs_m_xs = [point[1] for point in vs_m_pts]
vs_m_ys = [point[2] for point in vs_m_pts]
vs_m_interp = interpolate(vs_m_xs,vs_m_ys, BSplineOrder(3))

#VS/SkX boundary
vs_skx_pts_indices = [169,238,332,355] 
vs_skx_pts = [boundary_pts[106]]
for i in eachindex(vs_skx_pts_indices)
    push!(vs_skx_pts, boundary_pts[vs_skx_pts_indices[i]])
end
vs_skx_xs = [point[1] for point in vs_skx_pts]
vs_skx_ys = [point[2] for point in vs_skx_pts]
vs_skx_interp = interpolate(vs_skx_xs,vs_skx_ys, BSplineOrder(3))

#M upper boundary
m_pts_indices = [28,62,106]
m_pts = [boundary_pts[2]]
for i in eachindex(m_pts_indices)
    push!(m_pts, boundary_pts[m_pts_indices[i]])
end
m_xs = [point[1] for point in m_pts]
m_ys = [point[2] for point in m_pts]
m_interp = interpolate(m_xs,m_ys, BSplineOrder(3))

#F left boundary
l_pts_indices = [9,18] 
l_pts = [boundary_pts[2]]
for i in eachindex(l_pts_indices)
    push!(l_pts, boundary_pts[l_pts_indices[i]])
end
l_xs = [point[1] for point in l_pts]
l_ys = [point[2] for point in l_pts]
l_interp = interpolate(l_xs,l_ys, BSplineOrder(3))

#SkX lower left boundary
skxl_pts_indices = [18,28] 
skxl_pts = [(1.005,0.245)]
for i in eachindex(skxl_pts_indices)
    push!(skxl_pts, boundary_pts[skxl_pts_indices[i]])
end
skxl_xs = [point[1] for point in skxl_pts]
skxl_ys = [point[2] for point in skxl_pts]
skxl_interp = interpolate(skxl_xs,skxl_ys, BSplineOrder(3))

#SkX upper boundary
skxu_pts_indices = [50,119,180,240,355] 
skxu_pts = [(1.005,0.245)]
for i in eachindex(skxu_pts_indices)
    push!(skxu_pts, boundary_pts[skxu_pts_indices[i]])
end
skxu_xs = [point[1] for point in skxu_pts]
skxu_ys = [point[2] for point in skxu_pts]
skxu_interp = interpolate(skxu_xs,skxu_ys, BSplineOrder(4))

#2q, 2q' boundary
push!(sec_order_pts,skxu_pts[3])
so_xs = [point[1] for point in sec_order_pts]
so_ys = [point[2] for point in sec_order_pts]
so_interp = interpolate(so_xs,so_ys, BSplineOrder(3))

#FM  boundary
fm_pts_indices = [58,180] #corner pt added in next line
fm_pts = [(1.0,0.6)]
for i in eachindex(fm_pts_indices)
    push!(fm_pts, boundary_pts[fm_pts_indices[i]])
end
fm_xs = [point[1] for point in fm_pts]
fm_ys = [point[2] for point in fm_pts]
fm_interp = interpolate(fm_xs,fm_ys, BSplineOrder(3))

GLMakie.activate!() #Activate GLMakie or CairoMakie 
CairoMakie.activate!()
# Draw interpolants and scatter together
begin
    fig = Figure(size=(1920,1080),figure_padding=(10,50,20,30))
    ax = Axis(fig[1,1];
        xlabel = L"Δ",
        ylabel = L"h/J_1^-",
        xlabelsize = 72,
        ylabelsize = 72,
        xticklabelsize = 60,
        yticklabelsize = 60,
        xtickwidth = 2,
        ytickwidth = 2,
        xticksize = 10,
        yticksize = 10,
        spinewidth = 3,
        xticks = [0.95:0.05:1.35...],
        yticks = [0.0:0.1:0.6...],
        limits = (0.95,1.35,0.0,0.6),
        xgridvisible = false,
        ygridvisible = false
    )
    #Add color bands
    
    #CS band
    band!([0,1], [0, 0], [0.6, 0.6],color=(:lightgreen,0.9))

    #FP band
    xs = range(1,fm_xs[end],100)
    ys_u = [0.6 for _ in 1:100]
    band!(xs, fm_interp.(xs),ys_u, color = (:brown,0.5))
    xs = range(fm_xs[end],1.35,100)
    band!(xs, skxu_interp.(xs),ys_u, color = (:brown,0.5))

    #2q' -- 2 bands
    xs = range(fm_xs[1],so_xs[end],101)
    band!(xs, so_interp.(xs), fm_interp.(xs),color=(:yellow,0.5))
    xs = range(so_xs[end],fm_xs[end],101)
    band!(xs, skxu_interp.(xs), fm_interp.(xs),color=(:yellow,0.5))
    
    #2q -- 3 bands
    xs = range(1,1.005,51)
    band!(xs, l_interp.(xs), so_interp.(xs),color=(:orange,0.5))
    xs = range(skxl_xs[1],l_xs[end],51)
    band!(xs, l_interp.(xs), skxl_interp.(xs),color=(:orange,0.5))
    xs = range(skxu_xs[1],so_xs[end],101)
    band!(xs, skxu_interp.(xs), so_interp.(xs),color=(:orange,0.5))

    #SkX
    xs = range(skxu_xs[1],skxu_xs[end],150)
    ys_u = zeros(150)
    ys_l = zeros(150)
    for (i, x) in enumerate(xs)
        ys_u[i] = skxu_interp(x)
        if i<18 ys_l[i] = skxl_interp(x)
        elseif i<57 ys_l[i] = m_interp(x)
        else ys_l[i] = vs_skx_interp(x)
        end
    end 
    band!(xs, ys_l, ys_u,color=(:purple,0.3))
    
    #FL
    xs = range(l_xs[1],skxl_xs[end],80)
    ys_l = zeros(80)
    ys_u = zeros(80)
    for (i, x) in enumerate(xs)
        if i<41 ys_u[i] = l_interp(x)
        else ys_u[i] = skxl_interp(x)
        end
        ys_l[i] = m_interp(x)
    end 
    band!(xs, ys_l, ys_u,color=(:coral,0.5))

    #VS
    xs = range(vs_m_xs[1],vs_skx_xs[5],173)
    ys_u = zeros(173)
    ys_l = zeros(173)
    for (i, x) in enumerate(xs)
        if i<67 ys_u[i] = vs_m_interp(x)
        else ys_u[i] = vs_skx_interp(x)
        end
    end 
    band!(xs, ys_l, ys_u, color = (:lightcyan,0.9))

    #M
    xs = range(m_xs[1],m_xs[4],73)
    ys_u = zeros(73)
    ys_l = zeros(73)
    for (i, x) in enumerate(xs)
        ys_u[i] = m_interp(x)
        ys_l[i] = vs_m_interp(x)
    end 
    band!(xs, ys_l, ys_u, color = (:tomato,0.9))
    
    #Add phase boundary lines
    xs_1 = range(1.0, boundary_pts[106][1], 100)
    lines!(ax, xs_1, x -> vs_m_interp(x),color = :black, linewidth = 3)
    xs_2 = range(boundary_pts[106][1],1.35,100)
    lines!(ax, xs_2, x -> vs_skx_interp(x),color = :black, linewidth = 3)
    lines!(ax, xs_1, x -> m_interp(x),color = :black, linewidth = 3)
    xs_3 = range(1.0, boundary_pts[18][1],100)
    lines!(ax, xs_3, x -> l_interp(x), color = :black, linewidth = 3)
    xs_4 = range(1.005,boundary_pts[28][1],100)
    lines!(ax, xs_4, x -> skxl_interp(x), color = :black, linewidth = 3)
    xs_5 = range(1.005,1.35,100)
    lines!(ax, xs_5, x -> skxu_interp(x), color = :black, linewidth = 3)
    xs_6 = range(sec_order_pts[1][1],sec_order_pts[12][1],100)
    lines!(ax, xs_6, x -> so_interp(x), color = :black, linestyle=:dash, linewidth = 6)
    xs_7 = range(1.0,boundary_pts[180][1],100)
    lines!(ax, xs_7, x -> fm_interp(x), color = :black, linewidth = 3)
    vlines!(ax,[1.0], color = :black, linewidth = 3)

    text!(ax, 1.25, 0.1; text="VS", fontsize=72)
    text!(ax, 1.05, 0.09; text="M", fontsize=72)
    text!(ax, 1.25, 0.45; text="FP", fontsize=72)
    text!(ax, 1.15, 0.25; text="SkX", fontsize=72)
    text!(ax, 0.965, 0.3; text="CS", fontsize=72)
    text!(ax, 1.05, 0.462; text="2q'", fontsize=72)
    text!(ax, 1.03, 0.38; text="2q", fontsize=72)
    text!(ax, 1.018, 0.115; text="FL", fontsize=50)

    fig
end
save("low_energy_phase_diagram_v2.png",fig)