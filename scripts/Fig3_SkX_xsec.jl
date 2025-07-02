using DrWatson
@quickactivate "2025-Williams-Skyrmions_frustrated_quantum_dimers"
using GLMakie, Sunny, LinearAlgebra
using CairoMakie


#CP3 version of the proper screw spiral (YZ-S)
using CairoMakie
GLMakie.activate!()
CairoMakie.activate!()
    begin
    ps = [Point3f(x, 0, 0) for x in -9:1:3 for z in 0:1:1] #for y in -6:1:6
    ns = map(p -> 0.1 * Vec3f(p[2], p[3], p[1]), ps)

    for (i,point) in enumerate(ps)
        R = 5
        x,y,z = point[1], point[2], point[3]
        n = [0,0.5*cos(2*pi*x/12),0.5*sin(2*pi*x/12)]
        if z==1
            S = [-1*n[1]/sqrt(2), -1*n[2]/sqrt(2),0.25*(1+2*n[3])]
            
        else
            S = [1*n[1]/sqrt(2), 1*n[2]/sqrt(2),0.25*(1+2*n[3])]
        end
        ns[i] = S
        

    end
    ns[13] = [0,0,0]
    ns[14] = [0,0,0]
    #=ns[25] = [0,0,0]
    ns[26] = [0,0,0]=#
    lengths = norm.(ns)

    fig = Figure(size=(2400,2400))
    ax = LScene(fig[1,1];show_axis = false)
    Makie.Camera3D(ax.scene; projectiontype = Makie.Perspective,
                        eyeposition = Vec3(-9.5,-3,1.25),
                        lookat = Vec3(3,0,-1))
    arrows!(ax,
        ps, ns, fxaa=true, # turn on anti-aliasing
        color=lengths,
        linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.2),
        align = :center
    )
    #meshscatter!(ax, ps[1:2];color="black")#lengths[121:122])
    meshscatter!(ax, ps[13:14];color="black")
    fig

    #=scene = arrows(
        ps, ns, fxaa=true, # turn on anti-aliasing
        color=lengths,
        linewidth = 0.03, arrowsize = Vec3f(0.1, 0.12, 0.1),
        align = :center, axis=(type=Axis3,)
    )=#
end

#save("CP3spiral_full_v6.png",fig)

#CP3 version of the AFM-CS (at B = 0, can have any polarization plane, no dipole moment fluctuation)
#Currently giving YZ-S....
using CairoMakie
GLMakie.activate!()
CairoMakie.activate!()
    begin
    ps = [Point3f(x, 0, z) for x in -9:1:3 for z in 0:1:1] #for y in -6:1:6
    ns = map(p -> 0.1 * Vec3f(p[2], p[3], p[1]), ps)

    for (i,point) in enumerate(ps)
        R = 5
        x,y,z = point[1], point[2], point[3]
        n = [0,0.5*cos(2*pi*x/12),0.5*sin(2*pi*x/12)]
        if z==1
            S = [-1*n[1]/sqrt(2), -1*n[2]/sqrt(2),0.25*(1+2*n[3])]
            
        else
            S = [1*n[1]/sqrt(2), 1*n[2]/sqrt(2),0.25*(1+2*n[3])]
        end
        ns[i] = S
        

    end
    ns[13] = [0,0,0]
    ns[14] = [0,0,0]
    #=ns[25] = [0,0,0]
    ns[26] = [0,0,0]=#
    lengths = norm.(ns)

    fig = Figure(size=(2400,2400))
    ax = LScene(fig[1,1];show_axis = false)
    Makie.Camera3D(ax.scene; projectiontype = Makie.Perspective,
                        eyeposition = Vec3(-9.5,-3,1.25),
                        lookat = Vec3(3,0,-1))
    arrows!(ax,
        ps, ns, fxaa=true, # turn on anti-aliasing
        color=lengths,
        linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.2),
        align = :center
    )
    #meshscatter!(ax, ps[1:2];color="black")#lengths[121:122])
    meshscatter!(ax, ps[13:14];color="black")
    fig

    #=scene = arrows(
        ps, ns, fxaa=true, # turn on anti-aliasing
        color=lengths,
        linewidth = 0.03, arrowsize = Vec3f(0.1, 0.12, 0.1),
        align = :center, axis=(type=Axis3,)
    )=#
end



#Full CP3 skyrmion
using CairoMakie
GLMakie.activate!()
CairoMakie.activate!()
begin
    szs = Vector{Float64}(undef,11*11*2)
    ps = [Point3f(x, y, z) for x in -5:1:5 for y in -5:1:5 for z in 0:1:1]
    ns = map(p -> 0.1 * Vec3f(p[2], p[3], p[1]), ps)
    
    for (i,point) in enumerate(ps)
        R = 5
        x,y,z = point[1], point[2], point[3]
        if x > 0
            phi = atan(y/x)
        elseif x < 0
            phi = atan(y/x) + pi
        elseif y > 0
            phi = pi/2
        else
            phi = 3*pi/2
        end
        r = sqrt(x^2+y^2)
        if r > R
            theta = 0
        else
            theta = pi*(1-r/R)
        end
        #Switch b/t SkX1 and SkX2 by minus sign in front
        #n = [0.5*sin(theta)*cos(phi),0.5*sin(theta)*sin(phi),0.5*cos(theta)] #Neel type
        n = -[-0.5*sin(theta)*sin(phi),0.5*sin(theta)*cos(phi),0.5*cos(theta)] #Bloch type
        if z==1
            S = [-1*n[1]/sqrt(2), -1*n[2]/sqrt(2),0.25*(1+2*n[3])]
            
        else
            S = [1*n[1]/sqrt(2), 1*n[2]/sqrt(2),0.25*(1+2*n[3])]
        end
        ns[i] = S
        szs[i] = S[3]
        
    end
    #Comment out 121, 122 for SkX1
    #ns[121] = [0,0,0]
    #ns[122] = [0,0,0]
    #ns[61] = [0,0,0]
    lengths = norm.(ns)
    fig = Figure(size = (1000,1000))
    ax = LScene(fig[1,1];show_axis = false)
    Makie.Camera3D(ax.scene; projectiontype = Makie.Perspective,
                    eyeposition = Vec3(-1.4,-2,1.5),
                    lookat = Vec3(0,0,0))
    arrows!(ax,
        ps, ns, fxaa=true, # turn on anti-aliasing
        color=szs,
        linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.2),
        align = :center
    )
    for (i,point) in enumerate(ps)
        if norm(ns[i]) < 0.0001
            meshscatter!(ax, ps[i];color="black")#lengths[121:122])
        end
    end

    #Take out scatter for skx1
    #meshscatter!(ax, ps[121:122];color="black")#lengths[121:122])
    
    #Optional can add semi-transparent planes
    #poly!(Point3f[(-5, -5,0), (5, -5,0), (5, 5,0), (-5, 5,0)], color = RGBAf(0,0,1,0.2), strokecolor = :black, strokewidth = 1)
    #poly!(Point3f[(-5, -5,1), (5, -5,1), (5, 5,1), (-5, 5,1)], color = RGBAf(1,0,0,0.1), strokecolor = :black, strokewidth = 1)

    
    fig
end
#safesave("CP3skx1_fullv3.png",fig)



#########################################################################################
########## Fig 3a: SkX-I cross section
#########################################################################################

GLMakie.activate!()
CairoMakie.activate!()
#Cross section of CP3 skyrmion
begin
    ps = [Point3f(x, 0, z) for x in -5:1:5 for z in 0:1:1]
    ns = map(p -> 0.1 * Vec3f(p[2], p[3], p[1]), ps)
    szs = Vector{Float64}(undef,11*2)
    for (i,point) in enumerate(ps)
        R = 5
        x,y,z = point[1], point[2], point[3]
        if x > 0
            phi = atan(y/x)
        elseif x < 0
            phi = atan(y/x) + pi
        elseif y > 0
            phi = pi/2
        else
            phi = 3*pi/2
        end
        r = sqrt(x^2+y^2)
        if r > R
            theta = 0
        else
            theta = pi*(1-r/R)
        end
        #Switch b/t SkX1 and SkX2 by minus sign in front
        n = -[-0.5*sin(theta)*sin(phi),0.5*sin(theta)*cos(phi),0.5*cos(theta)] #Bloch type
        if z==1
            S = [-1*n[1]/sqrt(2), -1*n[2]/sqrt(2),0.25*(1+2*n[3])]
            
        else
            S = [1*n[1]/sqrt(2), 1*n[2]/sqrt(2),0.25*(1+2*n[3])]
        end
        ns[i] = S
        szs[i] = S[3]

    end
    
    #Use these for skx1
    #lengths = norm.(ns)
    ns[1]=[0,0,0]
    ns[1]=[0,0,0]
    ns[21]=[0,0,0]
    ns[22]=[0,0,0]
    

    fig = Figure(size=(1200,1100))
    ax = LScene(fig[1,1];show_axis = false)
    Makie.Camera3D(ax.scene; projectiontype = Makie.Perspective,
                #eyeposition = Vec3(-1.4,-2,1.5),
                #lookat = Vec3(0,0,0))
                eyeposition = Vec3(0,-5,1),#Vec3(-5.2,-3.49,1.25),#Vec3(-5.5,-3,1.25),
                lookat = Vec3(0,0,0))#Vec3(3,0,-1))
    arrows!(ax,
    ps, ns, fxaa=true, # turn on anti-aliasing
    color= szs,#lengths,
    linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.2),
    align = :center
    )
    #11:12 for skx2
    #1:2, 21:22 for skx1
    meshscatter!(ax, ps[1:2];color="black")
    meshscatter!(ax, ps[21:22];color="black")

    fig


end
save("CP3skx1_xsec_v5.pdf",fig)

#########################################################################################
########## Fig 3b: SkX-II cross section
#########################################################################################

GLMakie.activate!()
CairoMakie.activate!()
#Cross section of CP3 skyrmion
begin
    ps = [Point3f(x, 0, z) for x in -5:1:5 for z in 0:1:1]
    ns = map(p -> 0.1 * Vec3f(p[2], p[3], p[1]), ps)
    szs = Vector{Float64}(undef,11*2)
    for (i,point) in enumerate(ps)
        R = 5
        x,y,z = point[1], point[2], point[3]
        if x > 0
            phi = atan(y/x)
        elseif x < 0
            phi = atan(y/x) + pi
        elseif y > 0
            phi = pi/2
        else
            phi = 3*pi/2
        end
        r = sqrt(x^2+y^2)
        if r > R
            theta = 0
        else
            theta = pi*(1-r/R)
        end
        #Switch b/t SkX1 and SkX2 by minus sign in front
        n = [0.5*sin(theta)*cos(phi),0.5*sin(theta)*sin(phi),0.5*cos(theta)] #Neel type
        if z==1
            S = [-1*n[1]/sqrt(2), -1*n[2]/sqrt(2),0.25*(1+2*n[3])]
            
        else
            S = [1*n[1]/sqrt(2), 1*n[2]/sqrt(2),0.25*(1+2*n[3])]
        end
        ns[i] = S
        szs[i] = S[3]

    end
    
    #Use these for skx1
    #lengths = norm.(ns)
    #=ns[1]=[0,0,0]
    ns[1]=[0,0,0]
    ns[21]=[0,0,0]
    ns[22]=[0,0,0]
    =#

    fig = Figure(size=(1200,1100))
    ax = LScene(fig[1,1];show_axis = false)
    Makie.Camera3D(ax.scene; projectiontype = Makie.Perspective,
                #eyeposition = Vec3(-1.4,-2,1.5),
                #lookat = Vec3(0,0,0))
                eyeposition = Vec3(0,-5,1),#Vec3(-5.2,-3.49,1.25),#Vec3(-5.5,-3,1.25),
                lookat = Vec3(0,0,0))#Vec3(3,0,-1))
    arrows!(ax,
    ps, ns, fxaa=true, # turn on anti-aliasing
    color= szs,#lengths,
    linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.2),
    align = :center
    )
    #11:12 for skx2
    #1:2, 21:22 for skx1
    #meshscatter!(ax, ps[1:2];color="black")
    #meshscatter!(ax, ps[21:22];color="black")

    #poly!(Point3f[(-5, -2,-0.5), (5, -2,-0.5), (5, 2,-0.5), (-5, 2,-0.5)], color = RGBAf(0,0,1,0.2), strokecolor = :black, strokewidth = 1)
    #poly!(Point3f[(-5, -2,0.5), (5, -2,0.5), (5, 2,0.5), (-5, 2,0.5)], color = RGBAf(1,0,0,0.1), strokecolor = :black, strokewidth = 1)


    fig


    #=scene2 = arrows(
        ps, ns, fxaa=true, # turn on anti-aliasing
        color=lengths,
        linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.4),
        align = :center, axis=(type=Axis3,)
    )=#
end
save("CP3skx2_xsec_v5.pdf",fig)

