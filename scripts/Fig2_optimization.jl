using Sunny, LinearAlgebra, JLD2, SharedArrays
BLAS.set_num_threads(16)
Threads.nthreads()

#This initializes a dimer system with NN and NNN parallel and cross exchange interactions
#and applied magnetic field.
function entangled_model_from_alpha_b(; dims, J=1.0, Δ=1.2, α, B)
    latvecs = lattice_vectors(1,1,2,90,90,120)
    positions = [[0, 0, 0], [0, 0, 0.501]]
    cryst = Crystal(latvecs, positions; types = ["B", "T"])
    #units = Sunny.Units.theory
    sys = System(cryst, [1 => Moment(; s=1/2, g=-1), 2=>Moment(; s=1/2, g=-1)], :SUN; dims) #,units)

    J1 = -α*J
    J2 = 2α*J/(1 + √5)
    Jp1 = (Δ + 0.5)*J1
    Jc1 = (Δ - 0.5)*J1
    Jp2 = (Δ + 0.5)*J2
    Jc2 = (Δ - 0.5)*J2

    set_exchange!(sys, J, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, Jp1, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, Jp1, Bond(2, 2, [1, 0, 0]))
    set_exchange!(sys, Jp2, Bond(1, 1, [2, 1, 0]))
    set_exchange!(sys, Jp2, Bond(2, 2, [2, 1, 0]))
    set_exchange!(sys, Jc1, Bond(1, 2, [1, 0, 0]))
    set_exchange!(sys, Jc2, Bond(1, 2, [2, 1, 0]))
    set_field!(sys, (0,0,B))

    return Sunny.EntangledSystem(sys, [(1, 2)]), cryst
end

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
0.25*(sqrt(5)+1.0)/(sqrt(5)-1.0)

#From h_sat = 0.40, and 0.0<α<1.0, we have -0.73 < B < 1.0
#So choose 101 points for alpha and 347 points for B

begin
    N1 = 209
    N2 = 28
    Bs = range(0.0,1.0,N1)
    αs = range(0.6,1.0,N2)
    J=1.0
    Δ=1.2
    dims = (5,5,1)
    dims_full = (5,5,1,2)
    dims_ent = (5,5,1,1)
    energies = zeros(N1, N2)
    energies_buffer = zeros(N1, N2)
    coherent_buffers = [fill(zero(Sunny.CVec{4}),dims_ent) for _ in 1:N1, _ in 1:N2]
    syss_ent = [entangled_model_from_alpha_b(; dims, J, Δ, α=αs[j], B = (0,0,Bs[i])) for i in 1:N1, j in 1:N2]
end


#Initial set of randomize-minimize trials for each point in phase space
@time begin
    num_runs = 1 #number of randomize-minimize trials to conduct
    
    for trial in 1:num_runs
        Threads.@threads for idx in eachindex(energies) 
            sys_entangled = syss_ent[idx] 

            randomize_spins!(sys_entangled)
            optout = minimize_energy!(sys_entangled; maxiters=100_000)
            println("Optout: ", optout)

            energies_buffer[idx] = energy(sys_entangled)
                
            if energies[idx] > energies_buffer[idx] #store the new min energy and spin config 
                println("In trial $trial of initial rounds, we found a better energy $(energies_buffer[idx]) compared to $(energies[idx])")
                energies[idx] = energies_buffer[idx]
                coherent_buffers[idx] .= sys_entangled.sys.coherents
            else #reset the spin configuration to what it was before the most recent randomize-minimize iteration
                for site in Sunny.eachsite(sys_entangled)
                    set_coherent!(sys_entangled, coherent_buffers[idx][site], site)
                end
            end
        end
    end
end


# Perform the main iteration
begin
    total_runs = 2
    num_random_runs = 2 # Number of times to do a random initialization
    neighbor_radius = 2  # This will determine the number of trials from neighboring parameter values

    #Series of randomize-minimize trials but starting from the spin configuration of nearest neighbor with lowest energy
    for trial in 1:total_runs
        println("Trial $trial of $total_runs...")

        # First do specified number of random trials
        for _ in 1:num_random_runs
            Threads.@threads for idx in eachindex(energies) 
                sys_entangled = syss_ent[idx] 
                randomize_spins!(sys_entangled)
                optout = minimize_energy!(sys_entangled; maxiters=100_000)
                println("Optout: ", optout)

                energies_buffer[idx] = energy(sys_entangled)
                            
                if energies[idx] > energies_buffer[idx] #store the new min energy and spin config 
                    energies[idx] = energies_buffer[idx] # No need to copy. Single numbers are always copied by default.
                    coherent_buffers[idx] .= sys_entangled.sys.coherents
                else #reset the spin configuration to what it was before the most recent randomize-minimize iteration
                    for site in Sunny.eachsite(sys_entangled)
                        set_coherent!(sys_entangled, coherent_buffers[idx][site], site)
                    end
                end
            end
        end

        # Then do trials over neighboring parameter values 
        for j in 1:N2, i in 1:N1
            neighbors = neighbor_indices(i, j, N1, N2, neighbor_radius)
            for (x, y) in neighbors
                sys_entangled = syss_ent[i,j]

                # Set the initial spin configuration to the minimum found with a
                # neighboring set of parameter values.
                for site in Sunny.eachsite(sys_entangled)
                    set_coherent!(sys_entangled, coherent_buffers[x,y][site], site)
                end

                # Minimize the energy
                opt_out = minimize_energy!(sys_entangled; maxiters=100_000)
                (opt_out == -1) && println("Warning! Minimization failed.") # Good idea to check if the minimization worked
                energies_buffer[i,j] = energy(sys_entangled)
                
                if energies[i,j] > energies_buffer[i,j] #store the new min energy and spin config 
                    energies[i,j] = energies_buffer[i,j] # No need to copy. Single numbers are always copied by default.
                    coherent_buffers[i,j] .= sys_entangled.sys.coherents
                else #reset the spin configuration to what it was before the most recent randomize-minimize iteration
                    for site in Sunny.eachsite(sys_entangled)
                        set_coherent!(sys_entangled, coherent_buffers[i,j][site], site)
                    end
                end
            end
        end
    end
end

for j in 1:N2, i in 1:N1
    coherent_buffers[i,j] .= syss_ent[i,j].sys.coherents
end

save("entangled_coherents_run_alpha06_60.jld2","coherent_buffers",coherent_buffers)



####################################################################################################
