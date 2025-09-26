include("parameters.jl")
include("analysis.jl")
include("diagnostics.jl")
include("components/atmosphere.jl")
include("components/land.jl")
include("components/ocean.jl")
include("components/ice.jl")
include("coupled_simulation.jl")
include("monin_obukhov.jl")
include("postprocessing.jl")
import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities

export coupled_heat_equations, solve_coupler!, run_simulation


"""Sets u0 = u(t), t0 = t, and empties u."""
function reinit!(cs::Interfacer.CoupledSimulation, t)
    for sim in cs.model_sims
        Interfacer.reinit!(sim.integrator, sim.integrator.u, t0=t)
    end
end

function update_atmos_values!(cs)

end

function update_atmos_values!(cs, ::Val{:over_sea})
    bound_ocean_vals = get_field(cs.model_sims.ocean_sim, Val(:T_oce_sfc))
    ice_T = get_field(cs.model_sims.ice_sim, Val(:T_ice))
    update_field!(cs.model_sims.atmos_sim, bound_ocean_vals, ice_T)
    return bound_ocean_vals
end

function update_atmos_values!(cs, ::Val{:over_land})
    T_Ls = get_field(cs.model_sims.land_sim, Val(:T_Ls))
    update_field!(cs.model_sims.atmos_sim, T_Ls)
    return T_Ls
end

function update_ocean_values!(cs)
    bound_atmos_vals = get_field(cs.model_sims.atmos_sim, Val(:F_AO))
    update_field!(cs.model_sims.ocean_sim, bound_atmos_vals)
    return bound_atmos_vals
end

function update_ice_values!(cs)
    bound_atmos_vals = get_field(cs.model_sims.atmos_sim, Val(:T_atm_sfc))
    bound_ocean_vals = get_field(cs.model_sims.ocean_sim, Val(:T_oce_sfc))
    update_field!(cs.model_sims.ice_sim, bound_atmos_vals, bound_ocean_vals)
    return bound_atmos_vals, bound_ocean_vals
end

function update_land_values!(cs)
    bound_atmos_vals = get_field(cs.model_sims.atmos_sim, Val(:F_AL))
    update_field!(cs.model_sims.land_sim, bound_atmos_vals)
    return bound_atmos_vals
end

function set_time!(cs::Interfacer.CoupledSimulation, t)
    cs.dates.date[1] = TimeManager.current_date(cs, t)
end

function time_in_s(cs::Interfacer.CoupledSimulation)
    t = Dates.datetime2epochms(cs.dates.date[1])
    t0 = Dates.datetime2epochms(cs.dates.date0[1])
    return Int((t - t0) / 1e3)
end

function advance_simulation!(cs::Interfacer.CoupledSimulation, t_end::Float64, parallel::Bool)
    config = Val(cs.model_sims[1].params.config)
    if parallel
        FieldExchanger.step_model_sims!(cs.model_sims, t_end)
        update_atmos_values!(cs, config)
        update_ocean_values!(cs)
        bound_atmos_vals, bound_ocean_vals = update_ice_values!(cs)
        update_land_values!(cs)
    else
        Interfacer.step!(cs.model_sims.ice_sim, t_end)
        Interfacer.step!(cs.model_sims.ocean_sim, t_end)
        Interfacer.step!(cs.model_sims.land_sim, t_end)
        update_atmos_values!(cs, config)
        Interfacer.step!(cs.model_sims.atmos_sim, t_end)
        update_ocean_values!(cs)
        bound_atmos_vals, bound_ocean_vals = update_ice_values!(cs)
        update_land_values!(cs)
    end
    return bound_atmos_vals, bound_ocean_vals
end


function solve_coupler!(
    cs::Interfacer.CoupledSimulation;
    iterations::Int=1,
    parallel::Bool=false,
)
    (; Δt_cpl, tspan) = cs

    @info("Starting coupling loop")

    ϱ_A = nothing
    ϱ_O = nothing

    t0 = tspan[begin]

    for t = t0:Δt_cpl:(tspan[end]-Δt_cpl)
        @info("Current time: $t")
        set_time!(cs, t)
        reinit!(cs, t)

        iter = 1
        atmos_vals_list = []
        ocean_vals_list = []
        T_A_Γ = nothing
        T_O_Γ = nothing

        Checkpointer.checkpoint_sims(cs)

        for iter = 1:iterations
            @info("Current iter: $(iter)")
            if iter > 1
                Checkpointer.restart!(cs, cs.dirs.checkpoints, time_in_s(cs))
                reinit!(cs, t)
            end

            T_A_Γ_old = T_A_Γ
            T_O_Γ_old = T_O_Γ

            try
                T_A_Γ, T_O_Γ = advance_simulation!(cs, t + Δt_cpl, parallel)
            catch err
                if isa(err, UnstableError)
                    @warn("Unstable simulation!")
                    return NaN, NaN
                end
                rethrow()
            end

            push!(atmos_vals_list, T_A_Γ)
            push!(ocean_vals_list, T_O_Γ)

            if has_converged(T_A_Γ, T_A_Γ_old, T_O_Γ, T_O_Γ_old)
                @info "Termination criterion satisfied!"
            end
        end
        if iterations > 1
            ϱ_A = compute_ϱ_numerical(atmos_vals_list, parallel)
            ϱ_O = compute_ϱ_numerical(ocean_vals_list, parallel)
        end
    end
    set_time!(cs, tspan[end])
    return ϱ_A, ϱ_O
end

function run_simulation(
    physical_values;
    iterations=10,
    parallel=false,
)
    cs = get_coupled_sim(physical_values)
    ϱ_A, ϱ_O = solve_coupler!(
        cs,
        iterations=iterations,
        parallel=parallel,
    )
    return cs, ϱ_A, ϱ_O
end

function coupled_heat_equations(;
    iterations::Int=1,
    parallel::Bool=false,
    monin_obukhov::Bool=true,
    kwargs...,
)
    physical_values = SimulationParameters(; kwargs...)

    if monin_obukhov
        physical_values.C_H_AO = compute_C_H_AO(physical_values)
        physical_values.C_H_AI = compute_C_H_AI(physical_values)
        restore_physical_values!(physical_values)
    end

    cs, ϱ_A, ϱ_O = run_simulation(physical_values, iterations=iterations, parallel=parallel)
    return cs, ϱ_A, ϱ_O
end;
