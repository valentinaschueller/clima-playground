import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer

export SeaIce, thickness_rhs!, solve_surface_energy_balance, ice_init, get_field, update_field!

struct SeaIce{P,Y,D,I} <: Interfacer.SeaIceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

Interfacer.name(::SeaIce) = "SeaIce"

function thickness_rhs!(dh, h, cache, t)
    if cache.ice_model_type != :thickness_feedback
        dh.data = 0.0
        return
    end
    T_Is = solve_surface_energy_balance(cache)
    conduction = cache.k_I ./ h.data .* (T_Is .- cache.T_Ib)
    bottom_melt = cache.C_IO .* (cache.T_Ib .- cache.T_O)
    dh.data = (bottom_melt .- conduction) ./ (cache.ρ_I * cache.L)
end

function solve_surface_energy_balance(c)
    shortwave = (1 - c.alb_I) * c.SW_in
    longwave = c.ϵ * (c.LW_in - c.A)
    sensible_sfc = c.C_AI .* (c.T_A .- 273.15)
    conduction = (c.k_I / c.h_I_ini) * (c.T_Ib - 273.15)
    T_Is = (conduction + shortwave + longwave .+ sensible_sfc) ./ (c.k_I / c.h_I_ini + c.ϵ * c.B + c.C_AI)
    return min.(T_Is .+ 273.15, 273.15)
end

function ice_init(stepping, ics, space, cache)
    Δt = Float64(stepping.Δt_min) / stepping.nsteps_ice
    saveat = stepping.timerange[1]:stepping.Δt_min:stepping.timerange[end]

    ode_function = CTS.ClimaODEFunction((T_exp!)=thickness_rhs!)
    problem = SciMLBase.ODEProblem(ode_function, ics, stepping.timerange, cache)
    integrator = SciMLBase.init(
        problem,
        stepping.odesolver,
        dt=Δt,
        saveat=saveat,
        adaptive=false,
    )
    sim = SeaIce(cache, ics, space, integrator)
    return sim
end

Checkpointer.get_model_prog_state(sim::SeaIce) = sim.integrator.u

Interfacer.step!(sim::SeaIce, t) =
    Interfacer.step!(sim.integrator, t - sim.integrator.t)

Interfacer.reinit!(sim::SeaIce) = Interfacer.reinit!(sim.integrator)

function get_field(sim::SeaIce, ::Val{:T_ice})
    if sim.integrator.p.ice_model_type == :constant
        return sim.integrator.p.T_I_ini
    end
    return vec(solve_surface_energy_balance(sim.integrator.p))
end

function update_field!(sim::SeaIce, T_A, T_O)
    parent(sim.integrator.p.T_A) .= T_A
    parent(sim.integrator.p.T_O) .= T_O
end
