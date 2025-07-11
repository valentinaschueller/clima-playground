import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer
import ClimaDiagnostics as CD

export HeatEquationAtmos, heat_atm_rhs!, atmos_init, get_field, update_field!

struct HeatEquationAtmos{P,Y,D,I} <: Interfacer.AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

Interfacer.name(::HeatEquationAtmos) = "HeatEquationAtmos"

function heat_atm_rhs!(dT, T, p::SimulationParameters, t)
    F_sfc = p.a_I * p.C_AI * (T[1] - p.T_Is) + (1 - p.a_I) * flux_AO(T, p)

    # set boundary conditions
    C3 = CC.Geometry.WVector
    # note: F_sfc is converted to a Cartesian vector in direction 3 (vertical)
    bcs_bottom = CC.Operators.SetValue(C3(F_sfc))
    bcs_top = CC.Operators.SetValue(C3(Float64(0)))

    ## gradient and divergence operators needed for diffusion in tendency calc.
    ᶠgradᵥ = CC.Operators.GradientC2F()
    ᶜdivᵥ = CC.Operators.DivergenceF2C(bottom=bcs_bottom, top=bcs_top)

    @. dT.data = ᶜdivᵥ(p.k_A * ᶠgradᵥ(T.data)) / (p.ρ_A * p.c_A)
end


function atmos_init(odesolver, ics, space, p::SimulationParameters, output_dir)
    ode_function = CTS.ClimaODEFunction((T_exp!)=heat_atm_rhs!)
    problem = SciMLBase.ODEProblem(ode_function, ics, (p.t_0, p.t_max), p)

    Δt = p.Δt_min / p.n_t_A
    air_temperature = CD.DiagnosticVariable(;
        short_name="T_A",
        long_name="Air Temperature",
        standard_name="air_temperature",
        units="K",
        (compute!)=(out, Y, p, t) -> get_prognostic_data!(out, Y, p, t),
    )
    diagnostic_handler = CD.DiagnosticsHandler([get_diagnostic(air_temperature, space, p.Δt_min, output_dir)], ics, p, p.t_0, dt=Δt)
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

    sim = HeatEquationAtmos(p, ics, space, integrator)
    return sim
end

Checkpointer.get_model_prog_state(sim::HeatEquationAtmos) = sim.integrator.u

function Interfacer.step!(sim::HeatEquationAtmos, t)
    Interfacer.step!(sim.integrator, t - sim.integrator.t)
    check_stability(sim.integrator.u, sim.params.stable_range)
end

function flux_AO(T, p::SimulationParameters)
    return p.C_AO * (T[1] - p.T_O)
end

Interfacer.reinit!(sim::HeatEquationAtmos) = Interfacer.reinit!(sim.integrator)

function get_field(sim::HeatEquationAtmos, ::Val{:T_atm_sfc})
    return vec([fieldvec[1] for fieldvec in sim.integrator.sol.u])
end

function get_field(sim::HeatEquationAtmos, ::Val{:F_AO})
    return vec([flux_AO(fieldvec, sim.integrator.p) for fieldvec in sim.integrator.sol.u])
end

function Interfacer.add_coupler_fields!(coupler_field_names, ::HeatEquationAtmos)
    coupler_fields = [:F_AO, :T_atm_sfc]
    push!(coupler_field_names, coupler_fields...)
end

function update_field!(sim::HeatEquationAtmos, T_O, T_Is)
    sim.integrator.p.T_O = mean(T_O)
    sim.integrator.p.T_Is = mean(T_Is)
end
