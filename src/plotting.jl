using GLMakie

function dipole(Z, S)
    Sunny.Vec3(real(Z'*S*Z) for S in S)
end

function expanded_system(sys)
    x, y, z = sys.latsize
    @assert z == 1 "Can only expand a 2D system"
    cryst = Crystal(lattice_vectors(1.0, 1.0, 2.0, 90, 90, 120), [[0,0,0]], 1)
    units = Sunny.Units.theory
    sys_expanded = System(cryst, (x, y, 2), [SpinInfo(1; S=1/2, g=1)], :dipole; units)
    S = spin_matrices(1/2)
    S1, S2 = Sunny.to_product_space(S, S)

    Zs = sys.coherents[:, :, :]
    for j in 1:y, i in 1:x
        sys_expanded.dipoles[i, j, 1] = dipole(Zs[i, j, 1], S1)
        sys_expanded.dipoles[i, j, 2] = dipole(Zs[i, j, 1], S2)
    end

    return sys_expanded
end

function plot_expanded_spins(sys; kwargs...)
    sys_dip = expanded_system(sys)
    plot_spins(sys_dip; kwargs...)
end

function numbers_to_colors!(out::AbstractArray{Makie.RGBAf}, in::AbstractArray{<: Number}, colormap, colorrange)
    @assert length(out) == length(in)
    if isnothing(colorrange) || colorrange[1] >= colorrange[2] - 1e-8
        fill!(out, first(colormap))
    else
        cmin, cmax = colorrange
        len = length(colormap)
        map!(out, in) do c
            # If `cmin ≤ in[i] ≤ cmax` then `0.5 ≤ x ≤ len+0.5`
            x = (c - cmin) / (cmax - cmin) * len + 0.5
            # Round to integer and clip to range [1, len]
            colormap[max(min(round(Int, x), len), 1)]
        end
    end
    return nothing
end

set_alpha(c, alpha) = Makie.coloralpha(c, alpha)

function cell_center(ndims)
    if ndims == 3
        return [1, 1, 1] / 2
    elseif ndims == 2
        return [1, 1, 0] / 2
    else
        error("Unsupported `ndims=$ndims`.")
    end
end

function cell_diameter(latvecs, ndims)
    (a1, a2, a3) = eachcol(latvecs)
    if ndims == 3
        return max(norm(a1+a2+a3), norm(a1+a2-a3), norm(a1-a2+a3), norm(a1-a2-a3))
    elseif ndims == 2
        return max(norm(a1+a2), norm(a1-a2))
    else
        error("Unsupported `ndims=$ndims`.")
    end
end

function orient_camera!(ax; lookat, camshiftdir, upvector, camdist, orthographic)
    if orthographic
        eyeposition = lookat - camdist * camshiftdir
        projectiontype = Makie.Orthographic
    else
        eyeposition = lookat - 2.5 * camdist * camshiftdir
        projectiontype = Makie.Perspective
    end

    # Disable the key that would reset camera
    reset = false
    # Do not automatically "recenter" when adding objects
    center = false
    # No rotations on zoom
    zoom_shift_lookat = false
    # Mouse-drag rotations are SO(3) symmetric
    fixed_axis = false

    Makie.cam3d!(ax.scene; lookat, eyeposition, upvector, projectiontype, reset, center, fixed_axis,
                 zoom_shift_lookat, clipping_mode=:view_relative, near=0.01, far=100)
end

function orient_camera!(ax, latvecs; ghost_radius, ℓ0, orthographic, ndims, camdist=nothing)
    a1, a2, a3 = eachcol(latvecs)
    if ndims == 3
        lookat = (a1 + a2 + a3)/2
        camshiftdir = normalize(a1 + a2)
        upvector = normalize(a1 × a2)
    elseif ndims == 2
        lookat = (a1 + a2) / 2
        camshiftdir = -normalize(a1 × a2)
        upvector = normalize((a1 × a2) × a1)
    else
        error("Unsupported dimension: $ndims")
    end

    # The extra shift ℓ0 is approximately the nearest-neighbor distance
    camdist = @something camdist max(cell_diameter(latvecs, ndims)/2 + 0.8ℓ0, ghost_radius)

    orient_camera!(ax; lookat, camshiftdir, upvector, camdist, orthographic)
end

function all_offsets_within_distance(latvecs, rs, pt; min_dist=0, max_dist, nonzeropart=false)
    # box_lengths[i] represents the perpendicular distance between two parallel
    # boundary planes spanned by lattice vectors a_j and a_k (where indices j
    # and k differ from i)
    box_lengths = [a⋅b/norm(b) for (a,b) in zip(eachcol(latvecs), eachrow(inv(latvecs)))]
    n_max = round.(Int, max_dist ./ box_lengths, RoundUp)

    idxs = Int[]
    offsets = Vec3[]

    for (i, r) in enumerate(rs)
        for n1 in -n_max[1]:n_max[1], n2 in -n_max[2]:n_max[2], n3 in -n_max[3]:n_max[3]
            n = Vec3(n1, n2, n3)
            nonzeropart && iszero(n) && continue

            dist = norm(latvecs * (r + n - pt))
            if min_dist <= dist <= max_dist
                push!(idxs, i)
                push!(offsets, n)
            end
        end
    end

    return (idxs, offsets)
end

function characteristic_length_between_atoms(cryst::Crystal)
    # Detect if atom displacements are on a submanifold (aligned line or plane)
    ps = cryst.positions[1:end-1] .- Ref(cryst.positions[end])
    any_nonzero = map(1:3) do i
        any(p -> !iszero(p[i]), ps)
    end
    vecs = eachcol(cryst.latvecs)[findall(any_nonzero)]

    # Take nth root of appropriate hypervolume per atom
    if length(vecs) == 0
        ℓ = Inf                            # For a single atom, use ℓ0 below
    elseif length(vecs) == 1
        ℓ = norm(vecs[1]) / Sunny.natoms(cryst)  # Atoms aligned with single lattice vector
    elseif length(vecs) == 2
        ℓ = sqrt(norm(vecs[1] × vecs[2]) / Sunny.natoms(cryst))
    elseif length(vecs) == 3
        ℓ = cbrt(abs(det(cryst.latvecs)) / Sunny.natoms(cryst))
    else
        error("Internal error")
    end

    # An upper bound is the norm of the smallest lattice vector.
    ℓ0 = minimum(norm.(eachcol(cryst.latvecs)))

    return min(ℓ0, ℓ)
end

function cell_wireframe(latvecs, ndims)
    vecs = Makie.Point3f0.(eachcol(latvecs))
    ret = Tuple{Makie.Point3f0, Makie.Point3f0}[]

    origin = zero(Makie.Point3f0)

    if ndims == 3
        for j in 0:1, k in 0:1
            shift = j*vecs[2]+k*vecs[3]
            push!(ret, (origin+shift, vecs[1]+shift))
        end
        for i in 0:1, k in 0:1
            shift = i*vecs[1]+k*vecs[3]
            push!(ret, (origin+shift, vecs[2]+shift))
        end
        for i in 0:1, j in 0:1
            shift = i*vecs[1]+j*vecs[2]
            push!(ret, (origin+shift, vecs[3]+shift))
        end
    elseif ndims == 2
        for j in 0:1
            shift = j*vecs[2]
            push!(ret, (origin+shift, vecs[1]+shift))
        end
        for i in 0:1
            shift = i*vecs[1]
            push!(ret, (origin+shift, vecs[2]+shift))
        end
    end
    return ret
end


function plot_spins_alt!(ax, sys::System; notifier=Makie.Observable(nothing), arrowscale=1.0, stemcolor=:lightgray, color=:red,
    colorfn=nothing, colormap=:viridis, colorrange=nothing, show_cell=true, orthographic=false,
    ghost_radius=0, ndims=3, dims=nothing, camdist=nothing)
    isnothing(dims) || error("Use notation `ndims=$dims` instead of `dims=$dims`")

    sysdims = size(sys.dipoles)[1:3]
    if ndims == 2
        sysdims[3] == 1 || error("System not two-dimensional in (a₁, a₂)")
    elseif ndims == 1
        sysdims[[2,3]] == [1,1] || error("System not one-dimensional in (a₁)")
    end

    # Show bounding box of magnetic supercell in gray (this needs to come first
    # to set a scale for the scene in case there is only one atom).
    supervecs = sys.crystal.latvecs * diagm(Vec3(size(sys.dipoles)[1:3]))
    # Makie.linesegments!(ax, cell_wireframe(supervecs, ndims); color=:gray, linewidth=1.5)

    # Infer characteristic length scale between sites
    ℓ0 = characteristic_length_between_atoms(Sunny.orig_crystal(sys))

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
            (idxs, offsets) = all_offsets_within_distance(supervecs, rs, cell_center(ndims); max_dist=ghost_radius, nonzeropart=true)
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
        Makie.linesegments!(ax, cell_wireframe(Sunny.orig_crystal(sys).latvecs, ndims); color=:teal, linewidth=1.5)
        pos = [(3/4)*Makie.Point3f0(p) for p in eachcol(Sunny.orig_crystal(sys).latvecs)[1:ndims]]
        text = [Makie.rich("a", Makie.subscript(repr(i))) for i in 1:ndims]
        Makie.text!(ax, pos; text, color=:black, fontsize=20, font=:bold, glowwidth=4.0,
                    glowcolor=(:white, 0.6), align=(:center, :center), depth_shift=-1f0)
    end

    orient_camera!(ax, supervecs; ghost_radius, ℓ0, orthographic, ndims, camdist)

    return ax
end

function layer_from_bilayer(sys, layer=1)
    crystal = sys.sys.crystal
    na, nb, nc = size(sys.sys.coherents)[1:3]
    sys_layer = System(crystal, (na, nb, nc), [SpinInfo(1; S=1/2, g=2)], :dipole)
    for c in 1:nc, b in 1:nb, a in 1:na
        # set_dipole!(sys_layer, sys.sys_origin.dipoles[a, b, c, layer], (a, b, c, 1))
        sys_layer.dipoles[a, b, c, 1] = sys.sys_origin.dipoles[a, b, c, layer]
    end
    return sys_layer
end
