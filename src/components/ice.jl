import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer, Utilities
import ClimaDiagnostics as CD
import ClimaCore.MatrixFields: @name, FieldMatrixWithSolver, FieldMatrix

export SeaIce, thickness_rhs!, T_Is, ice_init, get_field, update_field!

struct SeaIce{P,Y,D,I} <: Interfacer.SeaIceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

function J_s(t, p, T_sfc)
    if isnothing(p.J_s)
        return p.C_AI * (p.T_A - T_sfc)
    end
    return p.J_s(t)
end

function F_O(t, p)
    if isnothing(p.F_O)
        return p.C_IO * (p.T_O - p.T_Ib)
    end
    return p.F_O(t)
end

function SW_net(t, p)
    return (1 - p.alb_I) * p.SW_in(t)
end

function LW_net(t, p, T_sfc)
    return p.ϵ * (p.LW_in(t) - (p.A + p.B * (T_sfc - 273)))
end

function thickness_rhs!(dh, h, p, t)
    if p.ice_model_type != :thickness_feedback
        dh.data = 0.0
        return
    end
    T_sfc = T_Is(p, h[1], t)
    F_A = SW_net(t, p) + LW_net(t, p, T_sfc) + J_s(t, p, T_sfc) + p.J_q(t)
    dh.data = -(F_A + F_O(t, p)) / (p.q_I)
end

function T_Is(p, h_I=nothing, t=0.0)
    if isnothing(h_I)
        h_I = p.h_I_ini
    end
    conduction = (p.k_I / h_I) * (p.T_Ib - 273)
    denominator = p.k_I / h_I + p.ϵ * p.B
    if isnothing(p.J_s)
        denominator += p.C_AI
    end
    T_eq = (conduction + SW_net(t, p) + LW_net(t, p, 273) + p.J_q(t) + J_s(t, p, 273)) / denominator
    return min(T_eq + 273.0, 273.0)
end

function Wfact(W, h, p, dtγ, t)
    jac = 0.0
    if p.ice_model_type == :thickness_feedback && (T_Is(p, h[1], t) < 273)
        c2 = (p.ϵ * p.B) / p.k_I
        if isnothing(p.J_s)
            c2 += p.C_AI / p.k_I
        end
        c1 = -(SW_net(t, p) + LW_net(t, p, p.T_Ib) + J_s(t, p, p.T_Ib) + p.J_q(t)) / (p.q_I)
        jac = -c1 * c2 / (1 + c2 * h[1])^2
    end
    @. W.matrix[@name(data), @name(data)] = dtγ * jac * (LinearAlgebra.I,) - (LinearAlgebra.I,)
end

function get_T_Is(out, h, p, t)
    T_sfc = T_Is(p, h[1], t)
    if isnothing(out)
        field = copy(h)
        parent(field) .= T_sfc
        return field.data
    end
    out .= T_sfc
end

function get_ice_odefunction(ics, ::Val{:implicit})
    jacobian = FieldMatrix((@name(data), @name(data)) => similar(ics.data, CC.MatrixFields.DiagonalMatrixRow{Float64}))
    T_imp! = SciMLBase.ODEFunction(thickness_rhs!; jac_prototype=FieldMatrixWithSolver(jacobian, ics), Wfact=Wfact)
    return CTS.ClimaODEFunction((T_imp!)=T_imp!)
end

function get_ice_odefunction(ics, ::Val{:explicit})
    return CTS.ClimaODEFunction((T_exp!)=thickness_rhs!)
end

function ice_init(odesolver, p, output_dir)
    context = Utilities.get_comms_context(Dict("device" => "auto"))
    space = CC.Spaces.PointSpace(context, CC.Geometry.ZPoint(0.0))
    field_h_I = CC.Fields.ones(space) .* p.h_I_ini
    ics = CC.Fields.FieldVector(data=field_h_I)

    ode_function = get_ice_odefunction(ics, Val(p.timestepping))
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
    sim = SeaIce(p, ics, space, integrator)
    return sim
end

Checkpointer.get_model_prog_state(sim::SeaIce) = sim.integrator.u

function Interfacer.step!(sim::SeaIce, t)
    Interfacer.step!(sim.integrator, t - sim.integrator.t)
end

function Interfacer.get_field(sim::SeaIce, ::Val{:T_ice})
    h_I = Interfacer.get_field(sim, Val(:h_I))
    if sim.integrator.p.ice_model_type == :constant
        return sim.integrator.p.T_I_ini .* ones(size(h_I))
    end
    return vec([T_Is(sim.integrator.p, h, sim.integrator.t) for h in h_I])
end

function Interfacer.get_field(sim::SeaIce, ::Val{:h_I})
    return vec([fieldvec[end] for fieldvec in sim.integrator.sol.u])
end

function Interfacer.add_coupler_fields!(coupler_field_names, ::SeaIce)
    coupler_fields = [:T_ice, :h_I]
    push!(coupler_field_names, coupler_fields...)
end

function Interfacer.update_field!(sim::SeaIce, T_A, T_O)
    sim.integrator.p.T_A = mean(T_A)
    sim.integrator.p.T_O = mean(T_O)
end
