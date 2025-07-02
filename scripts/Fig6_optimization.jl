
using DrWatson
@quickactivate "2025-Williams-Skyrmions_frustrated_quantum_dimers"
using Sunny, GLMakie, LinearAlgebra, JLD2

includet(srcdir("model.jl"))
includet(srcdir("helpers.jl"))


# Calculate the indices of points in a uniformly-spaced grid (dimension N1xN2)
# lying within a given radius of point (i, j).
function neighbor_indices(i, j, N1, N2, radius)
    neighbors = []
    for (x, y) in Iterators.product(1:N1, 1:N2)
        if !(i==x && j==y) && √((x-i)^2 + (y-j)^2) < radius
            push!(neighbors, (x, y))
        end
    end
    return neighbors
end

function triangular_effective_model_nnn(; dims=(20, 20, 1), J1=1.0, J2=0.1, h=0.0, Δ = 1.0)
    # Create lattice and system
    cryst = Crystal(lattice_vectors(1.0, 1.0, 2.0, 90, 90, 120), [[0,0,0]])
    #units = Sunny.Units.theory
    units = Units(:meV, :angstrom)
    sys = System(cryst, dims, [SpinInfo(1; S=1/2, g=1)], :SUN; units)

    ex1 = J1 * [1 0 0;
                0 1 0;
                0 0 Δ]
    ex2 = J2 * [1 0 0;
                0 1 0;
                0 0 Δ]
    set_exchange!(sys, ex1, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, ex2, Bond(1, 1, [1, 2, 0]))
    #set_external_field!(sys, [0, 0, h])
    set_field!(sys, [0, 0, h])
    return sys, cryst
end

begin
    N1 = 55 #determine dimension of array of phase space points to sample
    N2 = 40
    hs = range(0,0.55,N1)
    Δs = range(0.95,1.35,N2)
    dims = (5,5,1)
    energies = zeros(N1, N2)
    dipole_buffers = [fill(zero(Sunny.Vec3),dims) for _ in 1:N1, _ in 1:N2]
end


#Allocate and initalize the systems and corresponding energies
begin
    J1 = -1.0
    J2 = 2/(1+sqrt(5))
    #Initialize syss and allocate at same time
    syss = [triangular_effective_model_nnn(; dims, J1, J2, Δ = Δs[j], h = hs[i])[1] for i in 1:N1, j in 1:N2]
    for (i, sys) in enumerate(syss)
        dipole_buffers[i] .= sys.dipoles
        energies[i] = energy(sys)
    end

    energies_initial = copy(energies)
end

#Initial set of randomize-minimize trials for each point in phase space
begin
    num_runs = 1 #number of randomize-minimize trials to conduct
    
    for trial in 1:num_runs
        @time for j in 1:N2, i in 1:N1
        
            randomize_spins!(syss[i,j])
            # Note that minimize_energy! does not return the minimum energy. It returns
            # the number of iterations it took to converge, or -1 in convergence failed.
            opt_out = minimize_energy!(syss[i, j]; maxiters=100_000)
            (opt_out == -1) && println("Warning! Minimization failed.") # Good idea to check if the minimization worked


            energy_buffer = energy(syss[i,j])
                
             if energies[i,j] > energy_buffer #store the new min energy and spin config 
                energies[i,j] = energy_buffer
                dipole_buffers[i,j] .= syss[i,j].dipoles
                #spins[i,j] .= dipole_buffer[i,j] #needed??
            else #reset the spin configuration to what it was before the most recent randomize-minimize iteration
                for site in Sunny.eachsite(syss[i,j])
                    set_dipole!(syss[i,j], dipole_buffers[i,j][site],site)
                end
            end
        end
        
    end
end

# Perform the main iteration
begin
    total_runs = 5
    num_random_runs = 1#10 # Number of times to do a random initialization
    neighbor_radius = 4  # This will determine the number of trials from neighboring parameter values


    #min_h = 1
    #min_Δ = 1

    #Series of randomize-minimize trials but starting from the spin configuration of nearest neighbor with lowest energy
    for trial in 1:total_runs
        println("Trial $trial of $total_runs...")

        # First do specified number of random trials
        @time for _ in 1:num_random_runs
            for j in 1:N2, i in 1:N1
                randomize_spins!(syss[i,j])
                # Note that minimize_energy! does not return the minimum energy. It returns
                # the number of iterations it took to converge, or -1 in convergence failed.
                opt_out = minimize_energy!(syss[i,j]; maxiters=100_000)
                (opt_out == -1) && println("Warning! Minimization failed.") # Good idea to check if the minimization worked
            
                # This calculates the energy
                energy_buffer = energy(syss[i,j])
                            
                if energies[i,j] > energy_buffer #store the new min energy and spin config 
                    energies[i,j] = energy_buffer # No need to copy. Single numbers are always copied by default.
                    dipole_buffers[i,j] .= syss[i,j].dipoles
                else #reset the spin configuration to what it was before the most recent randomize-minimize iteration
                    for site in Sunny.eachsite(syss[i,j])
                        set_dipole!(syss[i,j], dipole_buffers[i,j][site],site)
                    end
                end
            end
        end

        # Then do trials over neighboring parameter values
        @time for j in 1:N2, i in 1:N1
            neighbors = neighbor_indices(i, j, N1, N2, neighbor_radius)
            for (x, y) in neighbors

                # Set the initial spin configuration to the minimum found with a
                # neighboring set of parameter values.
                for site in Sunny.eachsite(syss[i,j])
                    set_dipole!(syss[i,j], dipole_buffers[x,y][site], site)
                end

                # Minimize the energy
                opt_out = minimize_energy!(syss[i,j]; maxiters=100_000)
                (opt_out == -1) && println("Warning! Minimization failed.") # Good idea to check if the minimization worked
                energy_buffer = energy(syss[i,j])
                
                if energies[i,j] > energy_buffer #store the new min energy and spin config 
                    energies[i,j] = energy_buffer # No need to copy. Single numbers are always copied by default.
                    dipole_buffers[i,j] .= syss[i,j].dipoles
                else #reset the spin configuration to what it was before the most recent randomize-minimize iteration
                    for site in Sunny.eachsite(syss[i,j])
                        set_dipole!(syss[i,j], dipole_buffers[i,j][site],site)
                    end
                end
            end
        end
    end
end




# Take a look at the heatmap of energies before and after optimization
begin
    fig = Figure()
    ax1 = Axis(fig[1,1]; title="Original energies", xlabel="h", ylabel="Δ")
    ax2 = Axis(fig[1,2]; title="Optimized energies", xlabel="H", ylabel="Δ")

    cmin = min(minimum(energies_initial), minimum(energies))
    cmax = max(maximum(energies_initial), maximum(energies))
    colorrange = (cmin, cmax)
    heatmap!(ax1, hs, Δs, energies_initial; colorrange)
    heatmap!(ax2, hs, Δs, energies; colorrange)

    fig
end

for j in 1:N2, i in 1:N1
    dipole_buffers[i,j] .= syss[i,j].dipoles
end

#Save dipole configurations which can be analyzed in Fig6_low_energy_phase_diagram.jl
save("dipoles_from_2200_pt_PS.jld2","dipole_buffers",dipole_buffers)






