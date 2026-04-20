import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer
import ClimaDiagnostics as CD
import ClimaCore.MatrixFields: @name, ⋅, FieldMatrixWithSolver, FieldMatrix
import LinearAlgebra
import Statistics: mean

export HeatEquationIce, heat_ice_rhs!, init

struct HeatEquationIce{P,Y,D,I} <: Interfacer.SeaIceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

function Wfact_heat_ice(W, Y, p, dtγ, t)
    C3 = CC.Geometry.WVector
    ᶜdivᵥ = CC.Operators.DivergenceF2C()
    ᶠgradᵥ = CC.Operators.GradientC2F(;
        bottom=CC.Operators.SetValue(C3(0.0)),
        top=CC.Operators.SetValue(C3(0.0)),
    )
    div_matrix = CC.MatrixFields.operator_matrix(ᶜdivᵥ)
    grad_matrix = CC.MatrixFields.operator_matrix(ᶠgradᵥ)
    @. W.matrix[@name(data), @name(data)] = (dtγ * p.k_I / (p.ρ_I * p.c_I)) * div_matrix() ⋅ grad_matrix() - (LinearAlgebra.I,)
end

function heat_ice_rhs!(dT, T, p::SimulationParameters, t)
    FT = eltype(p)
    T_sfc = T_Is(p, p.h_I_ini, t)

    # boundary conditions
    bcs_top = CC.Operators.SetValue(T_sfc)
    bcs_bottom = CC.Operators.SetValue(p.T_Ib)

    ## gradient and divergence operators needed for diffusion in tendency calc.
    ᶠgradᵥ = CC.Operators.GradientC2F(top=bcs_top, bottom=bcs_bottom)
    ᶜdivᵥ = CC.Operators.DivergenceF2C()

    @. dT.data = ᶜdivᵥ(p.k_I * ᶠgradᵥ(T.data)) / (p.ρ_I * p.c_I)
end

function get_heat_odefunction(ics, ::Val{:implicit})
    FT = eltype(ics)
    jacobian = FieldMatrix((@name(data), @name(data)) => similar(ics.data, CC.MatrixFields.TridiagonalMatrixRow{FT}))
    T_imp! = SciMLBase.ODEFunction(heat_ice_rhs!; jac_prototype=FieldMatrixWithSolver(jacobian, ics), Wfact=Wfact_heat_ice)
    return CTS.ClimaODEFunction((T_exp!)=nothing, (T_imp!)=T_imp!)
end

function get_heat_odefunction(ics, ::Val{:explicit})
    return CTS.ClimaODEFunction((T_exp!)=heat_ice_rhs!)
end

function init(odesolver, p::SimulationParameters, output_dir, ::Val{:heat_ice})
    FT = eltype(p)
    space = get_vertical_space(FT(0), p.h_I_ini, p.n_I)
    field_ice = CC.Fields.ones(space) .* p.T_I_ini
    ics = CC.Fields.FieldVector(data=field_ice)

    ode_function = get_heat_odefunction(ics, Val(p.timestepping))
    problem = SciMLBase.ODEProblem(ode_function, ics, (p.t_0, p.t_0 + p.t_max), p)

    Δt = p.Δt_min / p.n_t_O
    sea_ice_temperature = CD.DiagnosticVariable(;
        short_name="T_I",
        long_name="Sea Ice Temperature",
        standard_name="sea_ice_temperature",
        units="K",
        (compute!)=(out, Y, p, t) -> get_prognostic_data!(out, Y, p, t),
    )
    diagnostic_handler = CD.DiagnosticsHandler(
        [get_diagnostic(sea_ice_temperature, space, p.Δt_min, output_dir)],
        ics,
        p,
        p.t_0,
        dt=Δt,
    )
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

    sim = HeatEquationIce(p, ics, space, integrator)
    return sim
end

function Interfacer.get_field(sim::HeatEquationIce, ::Val{:T_ice})
    return vec([fieldvec[end] for fieldvec in sim.integrator.sol.u])
end

function Interfacer.get_field(sim::HeatEquationIce, ::Val{:h_I})
    return vec([sim.integrator.p.h_I_ini])
end
