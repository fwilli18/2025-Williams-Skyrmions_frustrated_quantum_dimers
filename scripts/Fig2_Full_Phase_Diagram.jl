using DrWatson
@quickactivate "2025-Williams-Skyrmions_frustrated_quantum_dimers"
using GLMakie, BSplineKit, Sunny, LinearAlgebra, JLD2, Revise

includet(srcdir("model.jl"))
includet(srcdir("helpers.jl"))

begin
FM_bound = [(0.06,0.9425),(0.1,0.9025),(0.2,0.8075),(0.3,0.7175),(0.36,0.6675),(0.4,0.6325),(0.43,0.6125),(0.44,0.6225),(0.48,0.6775),(0.52,0.7325)]
fan_bound = [(0.43,0.6125),(0.44,0.6025),(0.46,0.5725),(0.48,0.5175), (0.50,0.43),(0.51,0.25),(0.52,0.164)]#(0.50,0.4025),(0.5175,0.28),(0.52,0.164)] #(0.52,0.24) is iffy
Skx2_bound = [(0.06,0.9275),(0.09,0.895),(0.1,0.8825),(0.2,0.7625),(0.22,0.7375),(0.24,0.7175),
                (0.3,0.6725),(0.36,0.6275),(0.4,0.6075),(0.44,0.5825),(0.46,0.5725)] #add in point at (0.46,???)(0.45,0.585)
Cb2_bound = [(0.09,0.895),(0.1,0.8775),(0.12,0.8475),(0.14,0.8175),(0.2,0.7175),(0.22,0.6775),(0.2225,0.6725),
                (0.24,0.6425),(0.28,0.5925),(0.32,0.5475),(0.36,0.5025),(0.38,0.4775)]#(0.24,0.6325),
                #(0.26,0.6025),(0.3,0.5525),(0.4,0.4225),(0.44,0.3675),(0.4625,0.34)] #(0.22,0.6575),
q1_bound2 = [(0.06,0.9075),(0.1,0.8425),(0.14,0.7825),(0.16,0.7525),(0.2,0.6925),(0.21,0.6775),(0.2225,0.6725)]
Skx1_bound = [(0.06,0.8975),(0.11,0.8125),(0.14,0.755),(0.2,0.6825),(0.21,0.6725)] #add in point at (0.22,??)
str2_bound = [(0.06,0.8925),(0.08,0.8525),(0.09,0.84),(0.1,0.8225),(0.11,0.8125)] #add in point at (0.12,??)
q1_bound1 = [(0.09,0.84),(0.1,0.8075),(0.12,0.7625),(0.14,0.7175),(0.2,0.5525),(0.22,0.48),(0.24,0.4025),
                (0.26,0.3075),(0.28,0.165),(0.30,0.0)]#(0.2875,0.0)] #add point in at (0.22,??) based on qpm point to match
qpm_bound = [(0.06,0.8875),(0.14,0.7075),(0.2,0.5475)] #add point in at (0.22,)
str1_bound = [(0.22,0.7375),(0.24,0.7125),(0.26,0.6825),(0.32,0.6075),(0.38,0.4775),(0.4,0.4425),(0.44,0.3775),
                (0.47,0.3275),(0.48,0.3075),(0.5,0.2675),(0.51,0.25)]#(0.5175,0.28)] #add point at beginning that links up to Skx2 at (0.24,??) and at end to (0.50,??)
#(0.225,0.7275)
str3_bound = [(0.32,0.6075),(0.34,0.5825),(0.36,0.5675),(0.38,0.5475),(0.40,0.5275),(0.42,0.5125),(0.44,0.4925),
                (0.46,0.4775),(0.48,0.4625),(0.50,0.43)]#(0.50,0.4025)]
Cb1_bound = [(0.2,0.6925),(0.21,0.6725),(0.22,0.6575),(0.26,0.6025),(0.3,0.5525),(0.34,0.4975),(0.38,0.4475),
                (0.4,0.4225),(0.44,0.3675),(0.46,0.3375),(0.47,0.3275)] #add point at beg and end


#Collect boundary points for interpolations
#First make the FM lower boundary
fm_xs = [point[1] for point in FM_bound]
fm_ys = [point[2] for point in FM_bound]
fm_interp = interpolate(fm_xs,fm_ys, BSplineOrder(4))

#FM-CS
fan_xs = [point[1] for point in fan_bound]
fan_ys = [point[2] for point in fan_bound]
fan_interp = interpolate(fan_xs,fan_ys, BSplineOrder(3))

#Skx2 lower bound
Skx2_xs = [point[1] for point in Skx2_bound]
Skx2_ys = [point[2] for point in Skx2_bound]
Skx2_interp = interpolate(Skx2_xs,Skx2_ys, BSplineOrder(6))

#Cb2 lower boundary
Cb2_xs = [point[1] for point in Cb2_bound]
Cb2_ys = [point[2] for point in Cb2_bound]
Cb2_interp = interpolate(Cb2_xs,Cb2_ys, BSplineOrder(3))

#q1_bound2
q1_2xs = [point[1] for point in q1_bound2]
q1_2ys = [point[2] for point in q1_bound2]
q1_2interp = interpolate(q1_2xs,q1_2ys, BSplineOrder(3))

#SkX1 lower boundary
Skx1_xs = [point[1] for point in Skx1_bound]
Skx1_ys = [point[2] for point in Skx1_bound]
Skx1_interp = interpolate(Skx1_xs,Skx1_ys, BSplineOrder(3))

#str2 lower bound
str2_xs = [point[1] for point in str2_bound]
str2_ys = [point[2] for point in str2_bound]
str2_interp = interpolate(str2_xs,str2_ys, BSplineOrder(2))

#q1_bound1
q1_1xs = [point[1] for point in q1_bound1]
q1_1ys = [point[2] for point in q1_bound1]
q1_1interp = interpolate(q1_1xs,q1_1ys, BSplineOrder(3))

#qpm
qpm_xs = [point[1] for point in qpm_bound]
qpm_ys = [point[2] for point in qpm_bound]
qpm_interp = interpolate(qpm_xs,qpm_ys, BSplineOrder(3))

#str1 lower bound
str1_xs = [point[1] for point in str1_bound]
str1_ys = [point[2] for point in str1_bound]
str1_interp = interpolate(str1_xs,str1_ys, BSplineOrder(3))

#str3 lower bound
str3_xs = [point[1] for point in str3_bound]
str3_ys = [point[2] for point in str3_bound]
str3_interp = interpolate(str3_xs,str3_ys, BSplineOrder(3))

#Cb1 lower bound
Cb1_xs = [point[1] for point in Cb1_bound]
Cb1_ys = [point[2] for point in Cb1_bound]
Cb1_interp = interpolate(Cb1_xs,Cb1_ys, BSplineOrder(3))
end

using CairoMakie
GLMakie.activate!()
CairoMakie.activate!()
# Draw interpolants and scatter together -- plot full phase diagram
begin
    fig = Figure(size = (1920,1920))#size = (1920,1017))
    ax = Axis(fig[1,1];
        xlabel = L"α",
        ylabel = L"B/|J|",
        xlabelsize = 50,
        ylabelsize = 50,
        xticklabelsize = 72,#36,
        yticklabelsize = 72,#36,
        #title = "Coupled Dimers Full Phase Diagram",
        titlesize = 42,
        limits = (0.06,0.52,0.0,1.0), #was up to 0.52
        xticks = [0.06,0.1,0.2,0.3,0.4,0.5],#(0.1:0.1:0.5),
        yticks = (0.0:0.1:1.0),
        ylabelrotation = 0,
        xgridvisible = false,
        ygridvisible = false
    )

    #=cmin = minimum(chi[3:(N1-3),:])
    cmax = 0.5*maximum(chi)
    colorrange = (cmin, cmax)
    hm = heatmap!(ax, αs,Bs, chi; colorrange)
    Colorbar(fig[:,end+1],hm, ticklabelsize = 24)=#
    #scatter!(ax, boundary_pts)

    #colsize!(fig.layout, 1, Aspect(1, 0.6183))
    #=scatter!(ax, FM_bound)
    scatter!(ax, fan_bound)
    scatter!(ax, Skx2_bound)
    scatter!(ax, Cb2_bound)
    scatter!(ax, q1_bound2)
    scatter!(ax, Skx1_bound)
    scatter!(ax, str2_bound)
    scatter!(ax, q1_bound1)
    scatter!(ax, qpm_bound)
    scatter!(ax, str1_bound)
    scatter!(ax, Cb1_bound)=#
    


    #Add color bands
    
    #QPM band bounded above by str2_bound (on [0.06,0.09]) and q1_bound1 on [0.09,0.30]
    ys_0 = zeros(97)
    xs_qpm = range(0.06,0.30,97)
    ys_qpm = zeros(97)
    for (i, x) in enumerate(xs_qpm)
        if i<14 ys_qpm[i] = str2_interp(x)
        else ys_qpm[i] = q1_1interp(x)
        end
    end
    #band!(xs_qpm,ys_0, ys_qpm, color = (:gray90,0.9))

    #2Q-III band
    xs_str = range(0.06,0.11,20)
    band!(xs_str, str2_interp.(xs_str), Skx1_interp.(xs_str),color=(:coral,0.9))

    #SkX-I band
    xs_skx1 = range(0.06,0.21,61)
    ys_u = zeros(61)
    for (i, x) in enumerate(xs_skx1)
        if i<58 ys_u[i] = q1_2interp(x)
        else ys_u[i] = Cb1_interp(x)
        end
    end 
    band!(xs_skx1, Skx1_interp.(xs_skx1), ys_u, color=(:lime,0.9))

    #YZ-S band
    xs = range(0.06,0.2225,66)
    ys_u = zeros(66)
    for (i, x) in enumerate(xs)
        if i<14 ys_u[i] = Skx2_interp(x)
        else ys_u[i] = Cb2_interp(x)
        end
    end 
    band!(xs, q1_2interp.(xs), ys_u, color=(:blue,0.3))
    
    #Skx-II band
    xs = range(0.06,0.46,121)
    ys_u = zeros(121)
    for (i, x) in enumerate(xs)
        if i<112 ys_u[i] = fm_interp(x)
        else ys_u[i] = fan_interp(x)
        end
    end 
    band!(xs, Skx2_interp.(xs), ys_u, color=(:yellow,0.9))

    #3Q-II
    xs = range(0.09,0.38,117)
    ys_u = zeros(117)
    for (i, x) in enumerate(xs)
        if i<53 ys_u[i] = Skx2_interp(x)
        else ys_u[i] = str1_interp(x)
        end
    end 
    band!(xs, Cb2_interp.(xs), ys_u,color=(:brown,0.3))

    #3Q-I
    xs = range(0.2,0.47,109)
    ys_u = zeros(109)
    for (i, x) in enumerate(xs)
        if i<10 ys_u[i] = q1_2interp(x)
        elseif i<73 ys_u[i] = Cb2_interp(x)
        else ys_u[i] = str1_interp(x)
        end
    end 
    band!(xs, Cb1_interp.(xs), ys_u,color=(:orange,0.5))

    #2Q-II
    xs = range(0.22,0.50,113)
    ys_u = zeros(113)
    ys_l = zeros(113)
    for (i, x) in enumerate(xs)
        if i<97 ys_u[i] = Skx2_interp(x)
        else ys_u[i] = fan_interp(x)
        end
        if i<41 ys_l[i] = str1_interp(x)
        else ys_l[i] = str3_interp(x)
        end
    end 
    band!(xs, ys_l, ys_u,color=(:purple,0.3))
    
    #2Q-I 
    xs = range(0.32,0.51,80)
    ys_u = zeros(80)
    for (i, x) in enumerate(xs)
        if i<76 ys_u[i] = str3_interp(x)
        else ys_u[i] = fan_interp(x)
        end
    end 
    band!(xs, str1_interp.(xs), ys_u,color=(:coral,0.5))

    #AFM CS
    xs = range(0.09,0.52,173)
    ys_u = zeros(173)
    ys_l = zeros(173)
    for (i, x) in enumerate(xs)
        if i<9 ys_u[i] = str2_interp(x)
        elseif i<49 ys_u[i] = Skx1_interp(x)
        elseif i<153 ys_u[i] = Cb1_interp(x)
        elseif i<170 ys_u[i] = str1_interp(x)
        else ys_u[i] = fan_interp(x)
        end
        if i<85 ys_l[i] = q1_1interp(x)
        end
    end 
    band!(xs, ys_l, ys_u, color = (:lightcyan,0.9))

    #FM-CS
    xs = range(0.43,0.52,73)
    ys_u = zeros(73)
    ys_l = zeros(73)
    for (i, x) in enumerate(xs)
        ys_u[i] = fm_interp(x)
        ys_l[i] = fan_interp(x)
    end 
    band!(xs, ys_l, ys_u, color = (:tomato,0.9))
    

    #arrows
    #=
    vlines!(ax, 0.08,ymin=0.725,ymax=0.855, color = :black)
    scatter!(0.08,0.85,marker=:utriangle, color =:black,markersize=15)
    linesegments!([0.1,0.15],[0.6,0.75],color=:black)
    scatter!(0.15,0.75,marker=:utriangle, color =:black,markersize=15, rotation = deg2rad(90))
    linesegments!([0.34,0.355],[0.445,0.49],color=:black)
    scatter!(0.355,0.49,marker=:utriangle, color =:black,markersize=15, rotation = deg2rad(45))
    =#

    xs_1 = range(0.06, 0.52, 100)
    lines!(ax, xs_1, x -> fm_interp(x), color = :black)
    xs_2 = range(0.43,0.52,100)
    lines!(ax, xs_2, x -> fan_interp(x), color = :black)
    xs_3 = range(0.06, 0.46,100)
    lines!(ax, xs_3, x -> Skx2_interp(x), color = :black)

    xs_4 = range(0.09,0.2225,50)#range(0.09,0.38,100)
    lines!(ax, xs_4, x -> Cb2_interp(x), color = :black, linestyle = :dash)
    xs_4b = range(0.2225,0.38,50)
    lines!(ax, xs_4b, x -> Cb2_interp(x), color = :black)

    xs_5 = range(0.06,0.2225,100)
    lines!(ax, xs_5, x -> q1_2interp(x), color = :black)
    xs_6 = range(0.06,0.21,100)
    lines!(ax, xs_6, x -> Skx1_interp(x), color = :black)

    xs_7 = range(0.06,0.09,50)#range(0.06,0.11,100)
    lines!(ax, xs_7, x -> str2_interp(x), color = :black, linestyle = :dash)
    xs_7b = range(0.09,0.11,50)#range(0.06,0.11,100)
    lines!(ax, xs_7b, x -> str2_interp(x), color = :black)

    xs_8 = range(0.09,0.3,100)
    lines!(ax, xs_8, x -> q1_1interp(x), color = :black, linestyle = :dash) #2nd order transition
    xs_9 = range(0.32,0.5,100)
    lines!(ax, xs_9, x -> str3_interp(x), color = :black)
    xs_10 = range(0.22,0.51,100)
    lines!(ax, xs_10, x -> str1_interp(x), color = :black)
    xs_11 = range(0.2,0.47,100)
    lines!(ax, xs_11, x -> Cb1_interp(x), color = :black)


    set_theme!(fonts = (; regular = "Times New Roman", bold = "Times New Roman Bold"))
    
    #text!(ax, 0.35, 0.8; text=rich("FP";font=:bold), fontsize=30)
    text!(ax, 0.35, 0.8; text="FP", font = :bold, fontsize=50)
    ####text!(ax, 0.48, 0.55; text="Fan", fontsize=22)
    text!(ax, 0.22, 0.74; text="SkX-II",font=:bold, fontsize=37)
    text!(ax, 0.17, 0.47; text="QPM",font=:bold, fontsize=50)
    text!(ax, 0.42, 0.525; text="2Q-II",font=:bold, fontsize=50)
    text!(ax, 0.42, 0.42; text="2Q-I",font=:bold, fontsize=50)
    text!(ax, 0.35, 0.2; text="AFM-CS", font=:bold, fontsize=50)
    text!(ax, 0.47, 0.57; text="FM-CS",font=:bold, fontsize=50)
    text!(ax, 0.27, 0.61; text="3Q-II",font=:bold, fontsize=42)
    text!(ax, 0.30, 0.40; text="3Q-I",font=:bold, fontsize=50)
    text!(ax, 0.12, 0.9; text="YZ-S",font=:bold, fontsize=50)
    text!(ax, 0.1, 0.56; text="SkX-I",font=:bold, fontsize=50)
    #text!(ax, 0.15, 0.74; text=L"SkXI", fontsize=14)
    text!(ax, 0.075, 0.7; text="2Q-III",font=:bold, fontsize=50)



    #xs = range(0.5, 1.0, 100)
    #lines!(ax, xs, x -> b2_interp(x))
    #lines!(ax, xs, x -> b3_interp(x))
    #arrows!((0.8,0.7),(0.0,0.15))

    #text!(ax, 0.25, 0.75; text="Phase I", fontsize=24)
    #text!(ax, 0.8, 0.55; text="Phase II", fontsize=24)
    #text!(ax, 0.5, 0.2; text="Phase III", fontsize=24)

    fig
end
#save("Full_Phase_Diagram_v6.png",fig)

