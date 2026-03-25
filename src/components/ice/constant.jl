import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer, Utilities
import ClimaDiagnostics as CD
import ClimaCore.MatrixFields: @name, FieldMatrixWithSolver, FieldMatrix
import Statistics: mean
import LinearAlgebra

export ConstantSeaIce, init

struct ConstantSeaIce{P,Y,D,I} <: Interfacer.SeaIceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end


function constant_rhs!(dh, h, p, t)
    dh.data = 0.0
end

function constant_Wfact(W, h, p, dtγ, t)
    jac = 0.0
    @. W.matrix[@name(data), @name(data)] = dtγ * jac * (LinearAlgebra.I,) - (LinearAlgebra.I,)
end

function init(odesolver, p::SimulationParameters, output_dir, ::Val{:constant})
    FT = eltype(p)
    context = Utilities.get_comms_context(Dict("device" => "auto"))
    space = CC.Spaces.PointSpace(context, CC.Geometry.ZPoint(FT(0.0)))
    field_h_I = CC.Fields.ones(space) .* p.h_I_ini
    ics = CC.Fields.FieldVector(data=field_h_I)

    ode_function = get_ice_odefunction(ics, constant_rhs!, constant_Wfact, Val(p.timestepping))
    problem = SciMLBase.ODEProblem(ode_function, ics, (p.t_0, p.t_max), p)
    Δt = p.Δt_min / p.n_t_I
    ice_thickness = CD.DiagnosticVariable(;
        short_name="h_I",
        long_name="Sea Ice Thickness",
        standard_name="sea_ice_thickness",
        units="m",
        (compute!)=(out, Y, p, t) -> get_prognostic_data!(out, Y, p, t),
    )
    ice_surface_temperature = CD.DiagnosticVariable(;
        short_name="T_Is",
        long_name="Sea Ice Surface Temperature",
        standard_name="sea_ice_surface_temperature",
        units="K",
        (compute!)=(out, Y, p, t) -> get_T_Is(out, Y, p, t),
    )
    diagnostics = [
        get_diagnostic(ice_thickness, space, p.Δt_min, output_dir),
        get_diagnostic(ice_surface_temperature, space, p.Δt_min, output_dir),
    ]
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
    return ConstantSeaIce(p, ics, space, integrator)
end

function Interfacer.get_field(sim::ConstantSeaIce, ::Val{:T_ice})
    h_I = Interfacer.get_field(sim, Val(:h_I))
    return sim.integrator.p.T_I_ini .* LinearAlgebra.ones(size(h_I))
end

function Interfacer.get_field(sim::ConstantSeaIce, ::Val{:h_I})
    return vec([sim.integrator.p.h_I_ini])
end


