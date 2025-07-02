function tl_dimer_model(; dims=(3, 1, 1), J=1.0, Jp1, Jp2, Jc1, Jc2)
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    positions = [[0, 0, 0], [0, 0, 0.501]]
    crystal = Crystal(latvecs, positions; types=["B", "T"])

    sys_origin = System(crystal, [1 => Moment(; s=1/2, g=-1), 2 => Moment(; s=1/2, g=-1)], :SUN; dims)

    set_exchange!(sys_origin, J, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys_origin, Jp1, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys_origin, Jp1, Bond(2, 2, [1, 0, 0]))
    set_exchange!(sys_origin, Jp2, Bond(1, 1, [2, 1, 0]))
    set_exchange!(sys_origin, Jp2, Bond(2, 2, [2, 1, 0]))
    set_exchange!(sys_origin, Jc1, Bond(1, 2, [1, 0, 0]))
    set_exchange!(sys_origin, Jc2, Bond(1, 2, [2, 1, 0]))

    dimers = [(1, 2)]
    sys = Sunny.EntangledSystem(sys_origin, dimers)

    return (; sys, crystal)
end

function tl_dimer_model_from_low_energy_params(; dims=(3, 1, 1), J=1.0, α, Δ) 
    (; Jp1, Jc1, Jp2, Jc2) = Js_from_α_Δ(; J, α, Δ)
    return tl_dimer_model(; dims, J, Jp1, Jp2, Jc1, Jc2) 
end

function Js_from_α_Δ(; J, α, Δ)
    J1 = -α*J
    J2 = 2α*J/(1 + √5)
    Jp1 = (Δ + 0.5)*J1
    Jc1 = (Δ - 0.5)*J1
    Jp2 = (Δ + 0.5)*J2
    Jc2 = (Δ - 0.5)*J2
    return (; Jp1, Jc1, Jp2, Jc2)
end