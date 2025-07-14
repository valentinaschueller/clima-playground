import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer
import ClimaDiagnostics as CD
import ClimaCore.MatrixFields: @name, FieldMatrixWithSolver, FieldMatrix

export SeaIce, thickness_rhs!, solve_surface_energy_balance, ice_init, get_field, update_field!

struct SeaIce{P,Y,D,I} <: Interfacer.SeaIceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

Interfacer.name(::SeaIce) = "SeaIce"

function thickness_rhs!(dh, h, p::SimulationParameters, t)
    if p.ice_model_type != :thickness_feedback
        dh.data = 0.0
        return
    end
    T_Is = solve_surface_energy_balance(p; h_I=vec(h.data))
    conduction = p.k_I ./ h .* (T_Is .- p.T_Ib)
    bottom_melt = p.C_IO * (p.T_Ib - p.T_O)
    @. dh = (bottom_melt - conduction) / (p.ρ_I * p.L)
end

function solve_surface_energy_balance(c; h_I=nothing)
    if isnothing(h_I)
        h_I = vec([c.h_I_ini])
    end
    T_Is = similar(h_I)
    shortwave = (1 - c.alb_I) * c.SW_in
    longwave = c.ϵ * (c.LW_in - c.A)
    sensible_sfc = c.C_AI * (c.T_A - 273.15)
    conduction = (c.k_I ./ h_I) .* (c.T_Ib .- 273.15)
    @. T_Is = (conduction + shortwave + longwave + sensible_sfc) / (c.k_I / h_I + c.ϵ * c.B + c.C_AI)
    return min.(T_Is .+ 273.15, 273.15)
end

function Wfact(W, h, p, dtγ, t)
    J = similar(h.data)
    if p.ice_model_type != :thickness_feedback
        @. J = 0.0
    elseif solve_surface_energy_balance(p; h_I=h)[1] < 273.15
        shortwave = (1 - p.alb_I) * p.SW_in
        longwave = p.ϵ * (p.LW_in - p.A - p.B * (p.T_Ib - 273.15))
        sensible = p.C_AI * (p.T_A - p.T_Ib)
        c1 = -(shortwave + longwave + sensible) / (p.ρ_I * p.L)
        c2 = (p.ϵ * p.B + p.C_AI) / p.k_I
        @. J = -c1 * c2 / (1 + c2 * h.data)^2
    else
        @info "Alternative Jacobian"
        c3 = -p.k_I / (p.ρ_I * p.L) * (273.15 - p.T_Ib)
        @. J = -c3 / h.data^2
    end
    @. W.matrix[@name(data), @name(data)] = dtγ * J * (LinearAlgebra.I,) - (LinearAlgebra.I,)
end

function get_T_Is(out, Y, p, t)
    field = copy(Y)
    field .= solve_surface_energy_balance(p; h_I=Y)
    if isnothing(out)
        return field.data
    else
        out .= field.data
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

function ice_init(odesolver, ics, space, p::SimulationParameters, output_dir)
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
    check_stability(get_field(sim, Val(:T_ice)), sim.params.stable_range)
end

Interfacer.reinit!(sim::SeaIce) = Interfacer.reinit!(sim.integrator)

function get_field(sim::SeaIce, ::Val{:T_ice})
    h_I = get_field(sim, Val(:h_I))
    if sim.integrator.p.ice_model_type == :constant
        return sim.integrator.p.T_I_ini .* ones(size(h_I))
    end
    return vec(solve_surface_energy_balance(sim.integrator.p; h_I=h_I))
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
