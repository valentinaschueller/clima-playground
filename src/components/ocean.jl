import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer
import ClimaDiagnostics as CD
import ClimaCore.MatrixFields: @name, ⋅, FieldMatrixWithSolver, FieldMatrix

export HeatEquationOcean, heat_oce_rhs!, ocean_init, get_field, update_field!

struct HeatEquationOcean{P,Y,D,I} <: Interfacer.OceanModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

function Wfact_oce(W, Y, p, dtγ, t)
    C3 = CC.Geometry.WVector
    ᶜdivᵥ = CC.Operators.DivergenceF2C(;
        bottom=CC.Operators.SetValue(C3(0.0)),
        top=CC.Operators.SetValue(C3(0.0)),
    )
    ᶠgradᵥ = CC.Operators.GradientC2F()
    div_matrix = CC.MatrixFields.operator_matrix(ᶜdivᵥ)
    grad_matrix = CC.MatrixFields.operator_matrix(ᶠgradᵥ)
    @. W.matrix[@name(data), @name(data)] = dtγ * div_matrix() ⋅ grad_matrix() - (LinearAlgebra.I,)
end

function heat_oce_rhs!(dT, T, p::SimulationParameters, t)
    F_sfc = p.a_I * p.C_IO * (p.T_Ib - T[end]) + (1 - p.a_I) * p.F_AO

    ## set boundary conditions
    C3 = CC.Geometry.WVector
    # note: F_sfc is converted to a Cartesian vector in direction 3 (vertical)
    bcs_top = CC.Operators.SetValue(C3(F_sfc))
    bcs_bottom = CC.Operators.SetValue(C3(Float64(0)))

    ## gradient and divergence operators needed for diffusion in tendency calc.
    ᶠgradᵥ = CC.Operators.GradientC2F()
    ᶜdivᵥ = CC.Operators.DivergenceF2C(bottom=bcs_bottom, top=bcs_top)

    @. dT.data = ᶜdivᵥ(p.k_O * ᶠgradᵥ(T.data)) / (p.ρ_O * p.c_O)
end

function get_oce_odefunction(ics, ::Val{:implicit})
    jacobian = FieldMatrix((@name(data), @name(data)) => similar(ics.data, CC.MatrixFields.TridiagonalMatrixRow{Float64}))
    T_imp! = SciMLBase.ODEFunction(heat_oce_rhs!; jac_prototype=FieldMatrixWithSolver(jacobian, ics), Wfact=Wfact_oce)
    return CTS.ClimaODEFunction((T_exp!)=nothing, (T_imp!)=T_imp!)
end

function get_oce_odefunction(ics, ::Val{:explicit})
    return CTS.ClimaODEFunction((T_exp!)=heat_oce_rhs!)
end

function ocean_init(odesolver, ics, space, p::SimulationParameters, output_dir)
    ode_function = get_oce_odefunction(ics, Val(p.timestepping))
    problem = SciMLBase.ODEProblem(ode_function, ics, (p.t_0, p.t_0 + p.t_max), p)

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

    saveat = p.t_0:p.Δt_min:p.t_max
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

function get_field(sim::HeatEquationOcean, ::Val{:T_oce_sfc})
    return vec([fieldvec[end] for fieldvec in sim.integrator.sol.u])
end

function update_field!(sim::HeatEquationOcean, F_AO)
    sim.integrator.p.F_AO = mean(F_AO)
end

function Interfacer.add_coupler_fields!(coupler_field_names, ::HeatEquationOcean)
    ocean_coupler_fields = [:T_oce_sfc,]
    push!(coupler_field_names, ocean_coupler_fields...)
end
