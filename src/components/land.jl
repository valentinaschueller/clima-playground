import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer

export Land, ∑sfc_flx!, solve_surface_energy_balance, land_init, get_field, update_field!

struct Land{P,Y,D,I} <: Interfacer.LandModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

Interfacer.name(::Land) = "Land"

function ∑sfc_flx!(dT, T, p::SimulationParameters, t)
    @. dT = p.F_AL
end

function land_init(odesolver, ics, space, p::SimulationParameters, output_dir)
    ode_function = CTS.ClimaODEFunction((T_exp!)=∑sfc_flx!)
    problem = SciMLBase.ODEProblem(ode_function, ics, (p.t_0, p.t_max), p)
    Δt = p.Δt_min / p.n_t_A
    sfc_temperature = CD.DiagnosticVariable(;
        short_name="T_Ls",
        long_name="Surface Temperature",
        standard_name="surface_temperature",
        units="K",
        (compute!)=(out, Y, p, t) -> get_prognostic_data!(out, Y, p, t),
    )
    diagnostics = [get_diagnostic(sfc_temperature, space, p.Δt_min, output_dir)]
    diagnostic_handler = CD.DiagnosticsHandler(diagnostics, ics, p, p.t_0, dt=Δt)
    diag_cb = CD.DiagnosticsCallback(diagnostic_handler)

    saveat = p.t_0:p.Δt_min:p.t_max

    integrator = SciMLBase.init(
        problem,
        odesolver,
        dt=Δt,
        saveat=saveat,
        adaptive=false,
        callback=SciMLBase.CallbackSet(diag_cb),
    )
    sim = Land(p, ics, space, integrator)
    return sim
end

Checkpointer.get_model_prog_state(sim::Land) = sim.integrator.u

function Interfacer.step!(sim::Land, t)
    Interfacer.step!(sim.integrator, t - sim.integrator.t)
    check_stability(get_field(sim, Val(:T_Ls)), sim.params.stable_range)
end

Interfacer.reinit!(sim::Land) = Interfacer.reinit!(sim.integrator)

function get_field(sim::Land, ::Val{:T_Ls})
    return vec([fieldvec[end] for fieldvec in sim.integrator.sol.u])
end

function Interfacer.add_coupler_fields!(coupler_field_names, ::Land)
    coupler_fields = [:T_Ls]
    push!(coupler_field_names, coupler_fields...)
end

function update_field!(sim::Land, flux_AL)
    sim.integrator.p.F_AL = mean(flux_AL)
end
