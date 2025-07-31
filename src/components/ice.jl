import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer
import ClimaDiagnostics as CD
import ClimaCore.MatrixFields: @name, FieldMatrixWithSolver, FieldMatrix

export SeaIce, thickness_rhs!, compute_T_Is, ice_init, get_field, update_field!

struct SeaIce{P,Y,D,I} <: Interfacer.SeaIceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

Interfacer.name(::SeaIce) = "SeaIce"

function thickness_rhs!(dh, h, p, t)
    t = (t / (3600 * 24 * 30)) % 12
    if p.ice_model_type != :thickness_feedback
        dh.data = 0.0
        return
    end
    T_Is = compute_T_Is(p, h[1], t)
    LW_out = p.ϵ * (p.A + p.B * (T_Is - 273))
    if isnothing(p.J_s)
        J_s = p.C_AI * (p.T_A - T_Is)
    else
        J_s = p.J_s(t)
    end
    if isnothing(p.F_O)
        F_O = p.C_IO * (p.T_O - p.T_Ib)
    else
        F_O = p.F_O(t)
    end
    F_A = (1 - p.alb_I) * p.SW_in(t) + p.LW_in(t) + J_s + p.J_q(t) - LW_out
    @. dh = -(F_A + F_O) / (p.q_I)
end

function compute_T_Is(p, h_I=nothing, t=0.0)
    if isnothing(h_I)
        h_I = p.h_I_ini
    end
    conduction = (p.k_I / h_I) .* (p.T_Ib - 273)
    if isnothing(p.J_s)
        T_eq = (conduction + (1 - p.alb_I) * p.SW_in(t) + p.ϵ * (p.LW_in(t) - p.A) + p.C_AI * (p.T_A - 273) + p.J_q(t)) / (p.k_I / h_I + p.ϵ * p.B + p.C_AI)
    else
        T_eq = (conduction + (1 - p.alb_I) * p.SW_in(t) + p.ϵ * (p.LW_in(t) - p.A) + p.J_s(t) + p.J_q(t)) / (p.k_I / h_I + p.ϵ * p.B)
    end
    return min(T_eq + 273.0, 273.0)
end

function Wfact(W, h, p, dtγ, t)
    jac = 0.0
    if p.ice_model_type == :thickness_feedback && (compute_T_Is(p, h[1], t) > 273)
        SW = (1 - p.alb_I) * p.SW_in(t)
        LW = p.ϵ * (p.LW_in(t) - p.A - p.B * (p.T_Ib - 273))
        if isnothing(p.J_s)
            J_s = p.C_AI * (p.T_A - p.T_Ib)
            c2 = (p.ϵ * p.B + p.C_AI) / p.k_I
        else
            J_s = p.J_s(t)
            c2 = (p.ϵ * p.B) / p.k_I
        end
        c1 = -(SW + LW + J_s + p.J_q(t)) / (p.q_I)
        jac = -c1 * c2 / (1 + c2 * h[1])^2
    end
    @. W.matrix[@name(data), @name(data)] = dtγ * jac * (LinearAlgebra.I,) - (LinearAlgebra.I,)
end

function get_T_Is(out, h, p, t)
    T_Is = compute_T_Is(p, h[1], t)
    if isnothing(out)
        field = copy(h)
        field .= T_Is
        return field.data
    else
        out .= T_Is
    end
end

function get_ice_odefunction(ics, ::Val{:implicit})
    jacobian = FieldMatrix((@name(data), @name(data)) => similar(ics.data, CC.MatrixFields.DiagonalMatrixRow{Float64}))
    T_imp! = SciMLBase.ODEFunction(thickness_rhs!; jac_prototype=FieldMatrixWithSolver(jacobian, ics), Wfact=Wfact)
    return CTS.ClimaODEFunction((T_imp!)=T_imp!)
end

function get_ice_odefunction(ics, ::Val{:explicit})
    return CTS.ClimaODEFunction((T_exp!)=thickness_rhs!)
end

function ice_init(odesolver, ics, space, p, output_dir)
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

Interfacer.reinit!(sim::SeaIce) = Interfacer.reinit!(sim.integrator)

function get_field(sim::SeaIce, ::Val{:T_ice})
    h_I = get_field(sim, Val(:h_I))
    if sim.integrator.p.ice_model_type == :constant
        return sim.integrator.p.T_I_ini .* ones(size(h_I))
    end
    return vec([compute_T_Is(sim.integrator.p, h, sim.integrator.t) for h in h_I])
end

function get_field(sim::SeaIce, ::Val{:h_I})
    return vec([fieldvec[end] for fieldvec in sim.integrator.sol.u])
end

function Interfacer.add_coupler_fields!(coupler_field_names, ::SeaIce)
    coupler_fields = [:T_ice, :h_I]
    push!(coupler_field_names, coupler_fields...)
end

function update_field!(sim::SeaIce, T_A, T_O)
    sim.integrator.p.T_A = mean(T_A)
    sim.integrator.p.T_O = mean(T_O)
end
