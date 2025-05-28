import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer

export SeaIce, thickness_rhs!, solve_surface_energy_balance, ice_init, get_field, update_field!

struct SeaIce{P,Y,D,I} <: Interfacer.SeaIceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

Interfacer.name(::SeaIce) = "SeaIce"

function thickness_rhs!(dh, h, cache, t)
    if cache.ice_model_type != :thickness_feedback
        dh.data = 0.0
        return
    end
    if cache.boundary_mapping == "cit"
        index = argmin(abs.(parent(CC.Fields.coordinate_field(cache.T_A)) .- t))
        T_O = parent(cache.T_O)[index]
    else
        T_O = vec(cache.T_O)
        index = nothing
    end
    T_Is = solve_surface_energy_balance(cache; h_I=vec(h.data), index=index)
    conduction = cache.k_I ./ h .* (T_Is .- cache.T_Ib)
    bottom_melt = cache.C_IO .* (cache.T_Ib .- T_O)
    @. dh = (bottom_melt - conduction) / (cache.ρ_I * cache.L)
end

function solve_surface_energy_balance(c; h_I=nothing, index=nothing)
    if isnothing(h_I)
        h_I = vec([c.h_I_ini])
    end
    if isnothing(index)
        T_A = vec(c.T_A)
    else
        T_A = parent(c.T_A)[index]
    end
    T_Is = zeros(size(h_I))
    shortwave = (1 - c.alb_I) * c.SW_in
    longwave = c.ϵ * (c.LW_in - c.A)
    sensible_sfc = c.C_AI .* (T_A .- 273.15)
    conduction = (c.k_I ./ h_I) .* (c.T_Ib .- 273.15)
    @. T_Is = (conduction + shortwave + longwave + sensible_sfc) / (c.k_I / h_I + c.ϵ * c.B + c.C_AI)
    return min.(T_Is .+ 273.15, 273.15)
end

function ice_init(odesolver, ics, space, p::SimulationParameters, output_dir)
    ode_function = CTS.ClimaODEFunction((T_exp!)=thickness_rhs!)
    problem = SciMLBase.ODEProblem(ode_function, ics, (p.t_0, p.t_0 + p.Δt_cpl), p)
    Δt = p.Δt_min / p.n_t_I
    air_temperature = CD.DiagnosticVariable(;
        short_name="h_I",
        long_name="Sea Ice Thickness",
        standard_name="sea_ice_thickness",
        units="m",
        (compute!)=(out, Y, p, t) -> get_prognostic_data!(out, Y, p, t),
    )
    diagnostic_handler = CD.DiagnosticsHandler([get_diagnostic(air_temperature, space, p.Δt_min, output_dir)], ics, p, p.t_0, dt=Δt)
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
    parent(sim.integrator.p.T_A) .= T_A
    parent(sim.integrator.p.T_O) .= T_O
end
