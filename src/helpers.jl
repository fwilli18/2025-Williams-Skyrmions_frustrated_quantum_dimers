using Sunny, LinearAlgebra, FFTW
using GLMakie
include(pkgdir(Sunny, "examples", "extra", "Plotting", "plotting2d.jl"))

include(pkgdir(Sunny, "ext", "PlottingExt.jl"))

includet(srcdir("model.jl"))

function su2_berry_curvature(n₁, n₂, n₃)
    res = n₁ ⋅ (n₂ × n₃) / (norm(n₁)*norm(n₂)*norm(n₃) + (n₁⋅n₂)*norm(n₃) + (n₁⋅n₃)*norm(n₂) + (n₂⋅n₃)*norm(n₁))
    (!isfinite(res)) && (res = 0.0)
    Ω = atan(res)
    return 2Ω
end

function su2_berry_curvature(n₁, n₂, n₃)
    res = n₁ ⋅ (n₂ × n₃) / (norm(n₁)*norm(n₂)*norm(n₃) + (n₁⋅n₂)*norm(n₃) + (n₁⋅n₃)*norm(n₂) + (n₂⋅n₃)*norm(n₁))
    (!isfinite(res)) && (res = 0.0)
    Ω = atan(res)
    return 2Ω
end

function plaquette_idcs(dims::Tuple{Int,Int,Int})
    dx, dy, dz = dims
    (dz != 1) && println("Warning: Ignoring lattice c vector.")
    Triple = Tuple{Int,Int,Int}
    idcs = Array{Tuple{Triple,Triple,Triple},4}(undef, 2, dx + 1, dy + 1, 1)
    for j ∈ 1:(dy+1)
        for i ∈ 1:(dx+1)
            idcs[1, i, j, 1] = (
                (mod1(i  , dx), mod1(j,   dy), 1),
                (mod1(i+1, dx), mod1(j+1, dy), 1),
                (mod1(i  , dx), mod1(j+1, dy), 1),
            )
            idcs[2, i, j, 1] = (
                (mod1(i,   dx), mod1(j,   dy), 1),
                (mod1(i+1, dx), mod1(j,   dy), 1),
                (mod1(i+1, dx), mod1(j+1, dy), 1),
            )
        end
    end
    return idcs
end

plaquette_idcs(x::Array{T,4}) where T = plaquette_idcs(size(x)[1:3])

function plaquette_map(f::Function, x::Array{T,4}) where T
    
    dims = size(x)
    @assert dims[3] == dims[4] == 1 "Multiple sites and multiple layers are not supported."
    dims = dims[1:3]
    out = zeros(Float64, 2, dims...)
    idcs_all = plaquette_idcs(x)
    for i in CartesianIndices(dims)
        idcs = idcs_all[1, i]
        out[1, i] = f(x[idcs[1]..., 1], x[idcs[2]..., 1], x[idcs[3]..., 1])
        idcs = idcs_all[2, i]
        out[2, i] = f(x[idcs[1]..., 1], x[idcs[2]..., 1], x[idcs[3]..., 1])
    end
    out
end

# Plotting functions
################################################################################
function plot_triangular_plaquettes!(ax, f, sys, n₁min, n₁max, n₂min, n₂max;
    colormap=(:PRGn, 0.5), colorrange=(-1.5, 1.5), texts=nothing, text_offset = (0.0, 0.0))


    v₁, v₂ = [1, 0, 0], [-1/2, √3/2, 0]
    #hidespines!(ax)
    #hidedecorations!(ax)

    # Plot panels
    plaq1(p) = Makie.Polygon(Point2f.([p, p+v₁+v₂, p+v₂]))
    plaq2(p) = Makie.Polygon(Point2f.([p, p+v₁, p+v₂+v₁]))

    v₀ = Point3f([0, 0, 0])

    χ = plaquette_map(f, sys.sys.coherents)
    pgons = Makie.Polygon[]
    color = Float64[]
    for r in n₁min:n₁max
        for c in n₂min:n₂max
            base = (r - 1) * v₁ + (c - 1) * v₂ + v₀
            push!(pgons, plaq1(base))
            push!(color, χ[1, r, c, 1, 1])
            push!(pgons, plaq2(base))
            push!(color, χ[2, r, c, 1, 1])
        end
    end
    println("maximum of Berry curvature ", maximum(abs.(color)))
    println("Total Berry curvature ", sum((color))/4pi)
    println("Total Berry curvature absolute value ", sum(abs.(color))/4pi)
    poly!(ax, pgons; color, colormap, colorrange)
    if !isnothing(texts)
        text!(ax, v₀[1] - text_offset[1], v₀[2] - text_offset[2]; text=texts[i], fontsize=14)
    end
end


function sun_berry_curvature(z₁, z₂, z₃)
    z₁, z₂, z₃ = normalize.((z₁, z₂, z₃))
    n₁ = z₁ ⋅ z₂
    n₂ = z₂ ⋅ z₃
    n₃ = z₃ ⋅ z₁
    return angle(n₁ * n₂ * n₃)
end




function plot_berry_curvature(esys::EntangledSystem; size=(600,200), kwargs...)
    plot_triangular_plaquettes(sun_berry_curvature, [esys.sys.coherents]; size, kwargs...)
end

function plot_berry_curvature(esys::Array{EntangledSystem}; size=(600, 200), kwargs...)
    plot_triangular_plaquettes(sun_berry_curvature, [esys.sys.coherents for esys in esys];
        size, kwargs...)
end

function structure_factor_2D(sys; alims=(-2/3, 2/3), blims=(-1/√3, 1/√3), c=0.0, npoints=200,contraction=:trace)
    # Make structure factor
    measure = ssf_perp(sys)
    ic = SampledCorrelationsStatic(sys,measure)#instant_correlations(sys)
    add_sample!(ic, sys)

    # Find wave vectors
    A = [1    0 0
        -1/2 1 0
        0    0 1]
    as = range(alims..., npoints)
    bs = range(blims..., npoints)
    qs_ortho = [[a, b, c] for a in as, b in bs]
    qs_rlu = [A * q for q in qs_ortho]

    # Retrieve S(q)
    formula = intensity_formula(ic, contraction)
    is = instant_intensities_interpolated(ic, qs_rlu, formula)

    return is, as, bs
end

function plot_static_structure_factor(frame; colormap=:Blues, colorrange1=(0, 100), colorrange2=(0, 100), savefig::Bool=false, savefig_path=nothing)
    # Reciprocal space base vector
    b₁, b₂ = 2π * [1, 1 / √3], 2π * [0, 2 / √3]
    # Linear system size

    L = size(frame)[1]
    # The resolution in momentum space
    δq₁, δq₂ = b₁ / L, b₂ / L

    # The basic parallelogram in momentum space
    paralg(p) = Makie.Polygon(Point2f.([p, p + δq₁, p + δq₁ + δq₂, p + δq₂]))

    # Get dipoles then apply fft
    cartesians = CartesianIndices((1:L, 1:L))
    Sα = zeros(ComplexF64, L, L, 3)
    for i in cartesians
        Sα[i, :] = frame[i, 1, 1]
    end
    fft!(Sα, [1, 2])
    
    # Get the structure factors for quadrupolars and octupolar
    Squad = abs2.(Sα[:, :, 1]) + abs2.(Sα[:, :, 2])
    Soctu = abs2.(Sα[:, :, 3])
    

    # Normalize to sites
    Squad = vcat(Squad...) ./ L^2
    Soctu = vcat(Soctu...) ./ L^2

    # Makes the data four times larger for later manipulations (see below)
    Squad = repeat(Squad, 4)
    Soctu = repeat(Soctu, 4)

    # Buffers for polygons to plot the structure factor
    # rt: right-top, lt: left-top, lb: left-bottom, and rb: right-bottom
    pgons_rt = Makie.Polygon[]
    pgons_lt = Makie.Polygon[]
    pgons_lb = Makie.Polygon[]
    pgons_rb = Makie.Polygon[]

    for i in cartesians
        base = Point2f((i[1] - 1) * δq₁ + (i[2] - 1) * δq₂)
        push!(pgons_rt, paralg(base))
        base = Point2f((i[1] - 1 - L) * δq₁ + (i[2] - 1) * δq₂)
        push!(pgons_lt, paralg(base))
        base = Point2f((i[1] - 1 - L) * δq₁ + (i[2] - 1 - L) * δq₂)
        push!(pgons_lb, paralg(base))
        base = Point2f((i[1] - 1) * δq₁ + (i[2] - 1 - L) * δq₂)
        push!(pgons_rb, paralg(base))
    end

    pgons = vcat(pgons_rt, pgons_lt, pgons_lb, pgons_rb)
    
    # The boundary coordinates of the Brillouin zone of the honeycomb lattice
    boundary_qx_coors = [4π / 3, 2π / 3, -2π / 3, -4π / 3, -2π / 3, 2π / 3, 4π / 3]
    boundary_qy_coors = [0, 2π / √3, 2π / √3, 0, -2π / √3, -2π / √3, 0]

    GLMakie.activate!()
    fig = Figure(backgroundcolor = :transparent)
    ax_l = Axis(fig[1, 1], aspect=DataAspect())
    hidespines!(ax_l)
    hidedecorations!(ax_l)
    poly!(ax_l, pgons; color=Squad, colormap, colorrange=colorrange1)
    lines!(ax_l, boundary_qx_coors, boundary_qy_coors, color=:black)
    text!(ax_l, 0, 2pi/sqrt(3)-0.7, text =L"\mathcal{S}^{\perp}(\textbf{q})",
         align = (:center, :center), fontsize = 40)
    GLMakie.xlims!(ax_l, -4π / 3, 4π / 3)
    GLMakie.ylims!(ax_l, -2π/√3, 2π/√3)
    ax_r  = Axis(fig[2, 1], aspect = DataAspect())
    hidespines!(ax_r)
    hidedecorations!(ax_r)
    poly!(ax_r, pgons; color=Soctu, colormap, colorrange=colorrange2)
    lines!(ax_r, boundary_qx_coors, boundary_qy_coors, color=:black)
    text!(ax_r, 0, 2pi/sqrt(3)-0.7, text =L"\mathcal{S}^{zz}(\textbf{q})",
         align = (:center, :center), fontsize = 40)
    GLMakie.xlims!(ax_r, -4π/3, 4π/3)
    GLMakie.ylims!(ax_r, -2π/√3, 2π/√3)
    fig

    if savefig
        CairoMakie.activate!()
        save(savefig_path, fig)
    end
    return fig
end



function magnetization(sys)
    M = Sunny.Vec3(0, 0, 0)
    for site in Sunny.eachsite(sys)
        M += magnetic_moment(sys, site)
    end
    return M / *(size(sys.dipoles)...)
end

function low_energy_params_from_full_params(; J, Jp1, Jc1, Jp2, Jc2, B)
    B = B isa Float64 ? B : B[3]

    h = B - J - 6*(Jp1 + Jp2 + Jc1 + Jc2)/4
    C = -(B + J) + 6*(Jp1 + Jp2 + Jc1 + Jc2)/8
    Δ1 = (Jp1 + Jc1)/(2*(Jp1 - Jc1))
    Δ2 = (Jp2 + Jc2)/(2*(Jp2 - Jc2))
    J1 = (Jp1 - Jc1)
    J2 = (Jp2 - Jc2)

    l = 2π/acos((-J1 -J2)/(2J2))
    println("Skyrmion length scale: ", l)

    return (; h, C, Δ1, Δ2, J1, J2)
end

function full_params_from_low_energy_params(; Δ, h, α, J=1.0)
    J1 = -α*J
    J2 = 2α*J/(1 + √5)
    Jp1 = (Δ + 0.5)*J1
    Jc1 = (Δ - 0.5)*J1
    Jp2 = (Δ + 0.5)*J2
    Jc2 = (Δ - 0.5)*J2
    B = h*J1 + J + 3*Δ*(J1+J2)
    return (; J, Jp1, Jc1, Jp2, Jc2, B)
end

#=function plot_spins!(ax, sys::System; notifier=Makie.Observable(nothing), arrowscale=1.0, stemcolor=:lightgray, color=:red,
    colorfn=nothing, colormap=:viridis, colorrange=nothing, show_cell=true, orthographic=false,
    ghost_radius=0, ndims=3, dims=nothing, show_sup = true)
isnothing(dims) || error("Use notation `ndims=$dims` instead of `dims=$dims`")
#warn_wglmakie()

if ndims == 2
sys.dims[3] == 1 || error("System not two-dimensional in (a₁, a₂)")
elseif ndims == 1
sys.dims[[2,3]] == [1,1] || error("System not one-dimensional in (a₁)")
end

# Show bounding box of magnetic supercell in gray (this needs to come first
# to set a scale for the scene in case there is only one atom).
if show_sup
    supervecs = sys.crystal.latvecs * diagm(Vec3(sys.dims))
    Makie.linesegments!(ax, cell_wireframe(supervecs, ndims); color=:gray, linewidth=1.5)
end
# Infer characteristic length scale between sites
ℓ0 = characteristic_length_between_atoms(orig_crystal(sys))

# Quantum spin-s, averaged over all sites. Will be used to normalize
# dipoles.
s0 = (sum(sys.Ns)/length(sys.Ns) - 1) / 2

# Parameters defining arrow shape
a0 = arrowscale * ℓ0
arrowsize = 0.4a0
linewidth = 0.12a0
lengthscale = 0.6a0
markersize = 0.8linewidth
arrow_fractional_shift = 0.6

# Positions in fractional coordinates of supercell vectors
rs = [supervecs \ global_position(sys, site) for site in eachsite(sys)]

for isghost in (false, true)
if isghost
alpha = 0.08
(idxs, offsets) = Sunny.all_offsets_within_distance(supervecs, rs, cell_center(ndims); max_dist=ghost_radius, nonzeropart=true)
else
alpha = 1.0
idxs = eachindex(rs)
offsets = [zero(Vec3) for _ in idxs]
end

# Every call to RGBf constructor allocates, so pre-calculate color
# arrays to speed animations
cmap_with_alpha = set_alpha.(Makie.to_colormap(colormap), Ref(alpha))
numeric_colors = zeros(size(sys.dipoles))
rgba_colors = zeros(Makie.RGBAf, size(sys.dipoles))

if isnothing(colorfn)
# In this case, we can precompute the fixed `rgba_colors` array
# according to `color`
if color isa AbstractArray
@assert length(color) == length(sys.dipoles)
if eltype(color) <: Number
   dyncolorrange = @something colorrange extrema(color)
   numbers_to_colors!(rgba_colors, color, cmap_with_alpha, dyncolorrange)
else
   map!(rgba_colors, color) do c
       set_alpha(Makie.to_color(c), alpha)
   end
end
else
c = set_alpha(Makie.to_color(color), alpha)
fill!(rgba_colors, c)
end
end

# These observables will be reanimated upon calling `notify(notifier)`.
vecs = Makie.Observable(Makie.Vec3f0[])
pts = Makie.Observable(Makie.Point3f0[])
pts_shifted = Makie.Observable(Makie.Point3f0[])
arrowcolor = Makie.Observable(Makie.RGBAf[])

Makie.on(notifier, update=true) do _
empty!.((vecs[], pts[], pts_shifted[], arrowcolor[]))

# Dynamically adapt `rgba_colors` according to `colorfn`
if !isnothing(colorfn)
numeric_colors .= colorfn.(CartesianIndices(sys.dipoles))
dyncolorrange = @something colorrange extrema(numeric_colors)
numbers_to_colors!(rgba_colors, numeric_colors, cmap_with_alpha, dyncolorrange)
end

for (site, n) in zip(idxs, offsets)
v = (lengthscale / s0) * vec(sys.dipoles[site])
pt = supervecs * (rs[site] + n)
pt_shifted = pt - arrow_fractional_shift * v
push!(vecs[], Makie.Vec3f0(v))
push!(pts[], Makie.Point3f0(pt))
push!(pts_shifted[], Makie.Point3f0(pt_shifted))
push!(arrowcolor[], rgba_colors[site])
end
# Trigger Makie redraw
notify.((vecs, pts, pts_shifted, arrowcolor))
# isnothing(color) || notify(arrowcolor)
end

# Draw arrows
linecolor = (stemcolor, alpha)
Makie.arrows!(ax, pts_shifted, vecs; arrowsize, linewidth, linecolor, arrowcolor, diffuse=1.15, transparency=isghost)

# Small sphere inside arrow to mark atom position
Makie.meshscatter!(ax, pts; markersize, color=linecolor, diffuse=1.15, transparency=isghost)
end

# Bounding box of original crystal unit cell in teal
if show_cell
Makie.linesegments!(ax, cell_wireframe(orig_crystal(sys).latvecs, ndims); color=:teal, linewidth=1.5)
pos = [(3/4)*Makie.Point3f0(p) for p in eachcol(orig_crystal(sys).latvecs)[1:ndims]]
text = [Makie.rich("a", Makie.subscript(repr(i))) for i in 1:ndims]
Makie.text!(ax, pos; text, color=:black, fontsize=20, font=:bold, glowwidth=4.0,
   glowcolor=(:white, 0.6), align=(:center, :center), depth_shift=-1f0)
end

orient_camera!(ax, supervecs; ghost_radius, ℓ0, orthographic, ndims)

return ax
end=#