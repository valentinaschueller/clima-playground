include("components/atmosphere.jl")
include("components/ocean.jl")
include("components/ice.jl")
include("parameters.jl")
include("coupled_simulation.jl")
include("convergence_factors.jl")
include("monin_obukhov.jl")
include("postprocessing.jl")
import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities

"""
When a file is saved, its always called the same thing. It has to be renamed for each iteration to not overwrite the old one.

**Arguments:**

-`cs::Interfacer.CoupledSimulation`: The coupled simulation.
-`iter::Int`: The current iteration.
-`time::Float64`: The current simulation time.

**Optional Keyword Arguments:**

-`reverse::Boolean`: The file name sometimes needs to be reverted back to the original name, default: `false`.

"""
function rename_files(cs::Interfacer.CoupledSimulation, iter, reverse=false)
    time = time_in_s(cs)
    for sim in cs.model_sims
        original_file = joinpath(
            cs.dirs.checkpoints,
            "checkpoint_" * Interfacer.name(sim) * "_$time.hdf5",
        )
        new_file = joinpath(
            cs.dirs.checkpoints,
            "checkpoint_" * Interfacer.name(sim) * "_$iter" * "_$time.hdf5",
        )
        if !reverse
            mv(original_file, new_file, force=true)
        else
            mv(new_file, original_file, force=true)
        end
    end
    pid = ClimaComms.mypid(cs.comms_ctx)
    original_file =
        joinpath(cs.dirs.checkpoints, "checkpoint_coupler_fields_$(pid)_$time.jld2")
    new_file = joinpath(
        cs.dirs.checkpoints,
        "checkpoint_coupler_fields_$(pid)" * "_$iter" * "_$time.jld2",
    )
    if !reverse
        mv(original_file, new_file, force=true)
    else
        mv(new_file, original_file, force=true)
    end
end

"""Resets integrator time."""
function reinit!(cs::Interfacer.CoupledSimulation, t)
    for sim in cs.model_sims
        Interfacer.reinit!(sim.integrator, sim.integrator.u, t0=t)
    end
end

"""Restarts simulations."""
function restart_sims!(cs::Interfacer.CoupledSimulation)
    @info "Reading checkpoint!"
    time = time_in_s(cs)
    rename_files(cs, 0, true)
    Checkpointer.restart!(cs, cs.dirs.checkpoints, time)
    rename_files(cs, 0)
end


function update_atmos_values!(cs, ice_T)
    ocean_states = copy(cs.model_sims.ocean_sim.integrator.sol.u)
    ocean_vals = extract_matrix(ocean_states, "oce")
    bound_ocean_vals = ocean_vals[end, :]
    ocean_T = mean(bound_ocean_vals)
    if cs.model_sims.atmos_sim.params.boundary_mapping == "mean"
        update_field!(cs.model_sims.atmos_sim, ocean_T, ice_T)
    else
        update_field!(cs.model_sims.atmos_sim, bound_ocean_vals, ice_T)
    end
    return ocean_vals, bound_ocean_vals
end

function update_ocean_values!(cs)
    atmos_states = copy(cs.model_sims.atmos_sim.integrator.sol.u)
    atmos_vals = extract_matrix(atmos_states, "atm")
    bound_atmos_vals = atmos_vals[1, :]
    if cs.model_sims.ocean_sim.params.boundary_mapping == "mean"
        atmos_T = mean(bound_atmos_vals)
        update_field!(cs.model_sims.ocean_sim, atmos_T)
    else
        update_field!(cs.model_sims.ocean_sim, bound_atmos_vals)
    end
    return atmos_vals, bound_atmos_vals
end

function set_time!(cs::Interfacer.CoupledSimulation, t)
    cs.dates.date[1] = TimeManager.current_date(cs, t)
end

function time_in_s(cs::Interfacer.CoupledSimulation)
    t = Dates.datetime2epochms(cs.dates.date[1])
    t0 = Dates.datetime2epochms(cs.dates.date0[1])
    return Int((t - t0) / 1e3)
end

"""
Runs the CoupledSimulation.

**Arguments:**

-`cs::Interfacer.CoupledSimulation`: A coupled simulation with atmosphere, ice and ocean.

**Optional Keyword Arguments:**

-`iterations::Int`: Number of iterations before the Schwarz iteration is terminated, default: `1`.
-`parallel::Boolean`: Whether to run the parallel or alternating Schwarz iteration, default: `false`.
"""
function solve_coupler!(
    cs::Interfacer.CoupledSimulation;
    iterations=1,
    parallel=false,
)
    (; Δt_cpl, tspan) = cs

    @info("Starting coupling loop")

    # Extract the initial value range for stability check
    lower_limit_temp, upper_limit_temp = initial_value_range(cs)

    ϱ_A = nothing
    ϱ_O = nothing

    for t = tspan[begin]:Δt_cpl:(tspan[end]-Δt_cpl)
        set_time!(cs, t)
        @info("Current time: $t")

        # Sets u0 = u(t), t0 = t, and empties u
        reinit!(cs, t)

        # Checkpoint to save initial values at this coupling step
        Checkpointer.checkpoint_sims(cs)
        rename_files(cs, 0)

        iter = 1
        atmos_vals_list = []
        ocean_vals_list = []
        bound_atmos_vals = nothing
        bound_ocean_vals = nothing

        ice_T = get_field(cs.model_sims.ice_sim, Val(:T_ice))

        for iter = 1:iterations
            @info("Current iter: $(iter)")
            if iter > 1
                restart_sims!(cs)
                reinit!(cs, t)
            end

            # Temperature values for the previous iteration.
            pre_bound_atmos_vals = bound_atmos_vals
            pre_bound_ocean_vals = bound_ocean_vals

            if parallel
                FieldExchanger.step_model_sims!(cs.model_sims, t + Δt_cpl)
                atmos_vals, bound_atmos_vals = update_atmos_values!(cs, ice_T)
                ocean_vals, bound_ocean_vals = update_ocean_values!(cs)
            else
                Interfacer.step!(cs.model_sims.ice_sim, t + Δt_cpl)
                Interfacer.step!(cs.model_sims.ocean_sim, t + Δt_cpl)
                ocean_vals, bound_ocean_vals = update_atmos_values!(cs, ice_T)

                Interfacer.step!(cs.model_sims.atmos_sim, t + Δt_cpl)
                atmos_vals, bound_atmos_vals = update_ocean_values!(cs)
            end

            atm_stable = is_stable(atmos_vals, upper_limit_temp, lower_limit_temp)
            oce_stable = is_stable(ocean_vals, upper_limit_temp, lower_limit_temp)
            if !(atm_stable && oce_stable)
                @info("Unstable simulation!")
                ϱ_A = atm_stable ? NaN : Inf
                ϱ_O = oce_stable ? NaN : Inf
                return ϱ_A, ϱ_O
            end

            push!(atmos_vals_list, bound_atmos_vals)
            push!(ocean_vals_list, bound_ocean_vals)

            Checkpointer.checkpoint_sims(cs)
            rename_files(cs, iter)

            if has_converged(
                bound_atmos_vals,
                pre_bound_atmos_vals,
                bound_ocean_vals,
                pre_bound_ocean_vals
            )
                break
            end
        end
        if iterations > 1
            ϱ_A, ϱ_O = compute_ϱ_numerical(atmos_vals_list, ocean_vals_list)
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

"""
Setup for running the coupled simulation and running it.

**Optional Keyword Arguments:**

-`iterations::Int`: Number of iterations before the Schwarz iteration is terminated, default: `1`.
-`parallel::Boolean`: Whether to run the parallel or alternating Schwarz iteration, default: `false`.
-`params::Dict`: Simulation parameters, default: `Dict{Symbol,Int}()`.

"""
function coupled_heat_equations(;
    iterations=1,
    parallel=false,
    monin_obukhov=true,
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
