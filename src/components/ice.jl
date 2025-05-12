import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer

struct SeaIce{P,Y,D,I} <: Interfacer.SeaIceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

Interfacer.name(::SeaIce) = "SeaIce"

function heat_ice_rhs!(dT, T, cache, t)
    # note: here we can add an update for the sea ice
    dT.data = 0
end

function solve_surface_energy_balance(c)
    k_I = 2.03
    α_I = 0.8
    h_I = 1.0
    A = 309.8
    B = 3.69
    ϵ = 0.98
    LW_in = 150.0
    SW_in = 200.0
    T_Ib = -1.8
    shortwave = (1 - α_I) * SW_in
    longwave = ϵ * (LW_in - A)
    sensible_sfc = c.C_AI * (c.T_A - 273.15)
    conduction = (k_I / h_I) * T_Ib
    T_sfc = (conduction + shortwave + longwave + sensible_sfc) / (k_I / h_I + ϵ * B + c.C_AI)
    return min(T_sfc + 273.15, 273.15)
end

function ice_temperature_feedback!(dt, T, cache, t)
    T_new = solve_surface_energy_balance(cache)
    dT.data = (T_new - T.data) / dt
end

function ice_init(stepping, ics, space, cache)
    Δt = Float64(stepping.Δt_min) / stepping.nsteps_ice
    saveat = stepping.timerange[1]:stepping.Δt_min:stepping.timerange[end]

    ode_function = CTS.ClimaODEFunction((T_exp!)=heat_ice_rhs!)
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

Interfacer.step!(sim::SeaIce, t) =
    Interfacer.step!(sim.integrator, t - sim.integrator.t)

Interfacer.reinit!(sim::SeaIce) = Interfacer.reinit!(sim.integrator)

get_field(sim::SeaIce, ::Val{:T_ice}) = sim.integrator.u[1]

function update_field!(sim::SeaIce, field)
    parent(sim.integrator.p.T_A)[1] = field
end
