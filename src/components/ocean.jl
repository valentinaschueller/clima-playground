import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer

export HeatEquationOcean, heat_oce_rhs!, ocean_init, get_field, update_field!

struct HeatEquationOcean{P,Y,D,I} <: Interfacer.OceanModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
Interfacer.name(::HeatEquationOcean) = "HeatEquationOcean"

function heat_oce_rhs!(dT, T, cache, t)
    F_sfc = (
        cache.a_I *
        cache.C_IO *
        (cache.T_Ib - T[end]) +
        (1 - cache.a_I) *
        cache.C_AO *
        (parent(cache.T_A)[1] - T[end])
    )

    ## set boundary conditions
    C3 = CC.Geometry.WVector
    # note: F_sfc is converted to a Cartesian vector in direction 3 (vertical)
    bcs_top = CC.Operators.SetValue(C3(F_sfc))
    bcs_bottom = CC.Operators.SetValue(C3(Float64(0)))

    ## gradient and divergence operators needed for diffusion in tendency calc.
    ᶠgradᵥ = CC.Operators.GradientC2F()
    ᶜdivᵥ = CC.Operators.DivergenceF2C(bottom=bcs_bottom, top=bcs_top)

    @. dT.data = ᶜdivᵥ(cache.k_O * ᶠgradᵥ(T.data)) / (cache.ρ_O * cache.c_O)
end

function ocean_init(odesolver, ics, space, p::SimulationParameters, output_dir)
    ode_function = CTS.ClimaODEFunction((T_exp!)=heat_oce_rhs!)
    problem = SciMLBase.ODEProblem(ode_function, ics, (p.t_0, p.t_0 + p.Δt_cpl), p)

    Δt = p.Δt_min / p.n_t_O
    sea_water_temperature = CD.DiagnosticVariable(;
        short_name="T_O",
        long_name="Sea Water Temperature",
        standard_name="sea_water_temperature",
        units="K",
        (compute!)=(out, Y, p, t) -> get_prognostic_data!(out, Y, p, t),
    )
    diagnostic_handler = CD.DiagnosticsHandler([get_diagnostic(sea_water_temperature, space, p.Δt_min, output_dir)], ics, p, p.t_0, dt=Δt)
    diag_cb = CD.DiagnosticsCallback(diagnostic_handler)

    saveat = p.t_0:p.Δt_min:p.Δt_cpl
    integrator = SciMLBase.init(
        problem,
        odesolver,
        dt=Δt,
        saveat=saveat,
        adaptive=false,
        callback=SciMLBase.CallbackSet(diag_cb),
    )

    sim = HeatEquationOcean(p, ics, space, integrator)
    return sim
end

Checkpointer.get_model_prog_state(sim::HeatEquationOcean) = sim.integrator.u

function Interfacer.step!(sim::HeatEquationOcean, t)
    Interfacer.step!(sim.integrator, t - sim.integrator.t)
    check_stability(sim.integrator.u, sim.params.stable_range)
end


Interfacer.reinit!(sim::HeatEquationOcean) = Interfacer.reinit!(sim.integrator)

function get_field(sim::HeatEquationOcean, ::Val{:T_oce_sfc})
    return vec([fieldvec[end] for fieldvec in sim.integrator.sol.u])
end

function update_field!(sim::HeatEquationOcean, T_A)
    parent(sim.integrator.p.T_A) .= vec([mean(T_A)])
end

function Interfacer.add_coupler_fields!(coupler_field_names, ::HeatEquationOcean)
    ocean_coupler_fields = [:T_oce_sfc,]
    push!(coupler_field_names, ocean_coupler_fields...)
end
