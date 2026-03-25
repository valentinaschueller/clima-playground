
import ClimaCoupler: Checkpointer, Interfacer, Utilities

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

function Interfacer.update_field!(sim::Interfacer.SeaIceModelSimulation, T_A, T_O)
    sim.integrator.p.T_A = mean(T_A)
    sim.integrator.p.T_O = mean(T_O)
end

function Interfacer.add_coupler_fields!(coupler_field_names, ::Interfacer.SeaIceModelSimulation)
    coupler_fields = [:T_ice, :h_I]
    push!(coupler_field_names, coupler_fields...)
end

function Interfacer.step!(sim::Interfacer.SeaIceModelSimulation, t)
    Interfacer.step!(sim.integrator, t - sim.integrator.t)
end

function T_Is(p, h_I, t)
    conduction = (p.k_I / h_I) * (p.T_Ib - 273)
    denominator = p.k_I / h_I + p.ϵ * p.B
    if isnothing(p.J_s)
        denominator += p.C_AI
    end
    T_eq = (conduction + SW_net(t, p) + LW_net(t, p, 273) + p.J_q(t) + J_s(t, p, 273)) / denominator
    return min(T_eq + 273.0, 273.0)
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

function get_ice_odefunction(ics, rhs!, wfact, ::Val{:implicit})
    FT = eltype(ics)
    jacobian = FieldMatrix((@name(data), @name(data)) => similar(ics.data, CC.MatrixFields.DiagonalMatrixRow{FT}))
    T_imp! = SciMLBase.ODEFunction(rhs!; jac_prototype=FieldMatrixWithSolver(jacobian, ics), Wfact=wfact)
    return CTS.ClimaODEFunction((T_imp!)=T_imp!)
end

function get_ice_odefunction(ics, rhs!, wfact, ::Val{:explicit})
    return CTS.ClimaODEFunction((T_exp!)=rhs!)
end

Checkpointer.get_model_prog_state(sim::Interfacer.SeaIceModelSimulation) = sim.integrator.u