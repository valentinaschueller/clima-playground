import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer
import ClimaDiagnostics as CD
import ClimaCore.MatrixFields: @name, ⋅, FieldMatrixWithSolver, FieldMatrix
import Statistics: mean

export HeatEquationAtmos, heat_atm_rhs!, atmos_init

struct HeatEquationAtmos{P,Y,D,I} <: Interfacer.AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

function Wfact_atm(W, Y, p, dtγ, t)
    # we can use homogeneous BC here, only the type is relevant to determine W
    ᶠgradᵥ = CC.Operators.GradientC2F(
        bottom=CC.Operators.SetValue(0.0),
        top=CC.Operators.SetGradient(CC.Geometry.WVector(0.0)),
    )
    ᶜdivᵥ = CC.Operators.DivergenceF2C()
    div_matrix = CC.MatrixFields.operator_matrix(ᶜdivᵥ)
    grad_matrix = CC.MatrixFields.operator_matrix(ᶠgradᵥ)
    @. W.matrix[@name(data), @name(data)] = (dtγ * p.k_A / (p.ρ_A * p.c_A)) * div_matrix() ⋅ grad_matrix() - (LinearAlgebra.I,)
    return nothing
end

function heat_atm_rhs!(dT, T, p::SimulationParameters, t)
    FT = eltype(p)

    # Dirichlet BC
    bcs_bottom = CC.Operators.SetValue(p.T_O)
    # Neumann BC
    bcs_top = CC.Operators.SetGradient(
        CC.Geometry.WVector(FT(0))
    )

    ## gradient and divergence operators needed for diffusion in tendency calc.
    ᶠgradᵥ = CC.Operators.GradientC2F(bottom=bcs_bottom, top=bcs_top)
    ᶜdivᵥ = CC.Operators.DivergenceF2C()

    @. dT.data = ᶜdivᵥ(p.k_A * ᶠgradᵥ(T.data)) / (p.ρ_A * p.c_A)
end

function get_atm_odefunction(ics, ::Val{:implicit})
    FT = eltype(ics)
    jacobian = FieldMatrix((@name(data), @name(data)) => similar(ics.data, CC.MatrixFields.TridiagonalMatrixRow{FT}))
    T_imp! = SciMLBase.ODEFunction(heat_atm_rhs!; jac_prototype=FieldMatrixWithSolver(jacobian, ics), Wfact=Wfact_atm)
    return CTS.ClimaODEFunction((T_exp!)=nothing, (T_imp!)=T_imp!)
end

function get_atm_odefunction(ics, ::Val{:explicit})
    return CTS.ClimaODEFunction((T_exp!)=heat_atm_rhs!)
end


function atmos_init(odesolver, p::SimulationParameters, output_dir)
    FT = eltype(p)
    space = get_vertical_space(FT(0.0), p.h_A, p.n_A)
    field_atm = CC.Fields.ones(space) .* p.T_A_ini
    ics = CC.Fields.FieldVector(data=field_atm)

    ode_function = get_atm_odefunction(ics, Val(p.timestepping))
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

function flux_AI(T, p::SimulationParameters)
    return p.C_AI * (T[1] - p.T_Is)
end

function flux_AO(T, p::SimulationParameters)
    return p.C_AO * (T[1] - p.T_O)
end

function surface_gradient(T, p::SimulationParameters)
    ᶠgradᵥ = CC.Operators.GradientC2F(;
        bottom=CC.Operators.SetValue(p.T_O),
        top=CC.Operators.SetGradient(CC.Geometry.WVector(0.0)),
    )
    # cast gradient to WVector to apply multiplication with Δz
    gradient_values = p.k_A.* CC.Geometry.WVector.(ᶠgradᵥ.(T.data))

    surface_gradient = parent(gradient_values)[1]
    return surface_gradient
end

function Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:T_atm_sfc})
    return vec([fieldvec[1] for fieldvec in sim.integrator.sol.u])
end

function Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:F_AO})
    flux_vector = vec([surface_gradient(fieldvec, sim.integrator.p) for fieldvec in sim.integrator.sol.u])
    return flux_vector
end

function Interfacer.add_coupler_fields!(coupler_field_names, ::HeatEquationAtmos)
    coupler_fields = [:F_AO, :T_atm_sfc]
    push!(coupler_field_names, coupler_fields...)
end

function Interfacer.update_field!(sim::HeatEquationAtmos, T_O, T_Is)
    sim.integrator.p.T_O = mean(T_O)
    sim.integrator.p.T_Is = mean(T_Is)
end
