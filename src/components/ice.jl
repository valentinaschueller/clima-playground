import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer

struct ConstantIce{P,Y,D,I} <: Interfacer.SeaIceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

Interfacer.name(::ConstantIce) = "ConstantIce"

function heat_ice_rhs!(dT, T, cache, t) # Something that has to be implemented right?
    dT.ice = 0
end

function ice_init(stepping, ics, space, cache)
    Δt = Float64(stepping.Δt_min) / stepping.nsteps_ice

    ode_function = CTS.ClimaODEFunction((T_exp!)=heat_ice_rhs!)
    problem = SciMLBase.ODEProblem(ode_function, ics, stepping.timerange, cache)
    integrator = SciMLBase.init(
        problem,
        stepping.odesolver,
        dt=Δt,
        saveat=Float64(stepping.Δt_coupler),
        adaptive=false,
    )
    sim = ConstantIce(cache, ics, space, integrator)
    return sim
end

Interfacer.step!(sim::ConstantIce, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t)
Interfacer.reinit!(sim::ConstantIce) = Interfacer.reinit!(sim.integrator)

get_field(sim::ConstantIce, ::Val{:T_ice}) = sim.integrator.u[1]
function update_field!(sim::ConstantIce, ::Val{:T_ice}, field)
    sim.params.T_ice = field
end