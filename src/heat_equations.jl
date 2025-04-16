include("components/atmosphere.jl")
include("components/ocean.jl")
include("components/ice.jl")
include("plotting.jl")
include("coupled_simulation.jl")
include("convergence_factors.jl")
include("postprocessing.jl")
include("parameters.jl")
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
function rename_files(cs::Interfacer.CoupledSimulation, iter, time, reverse=false)
    for sim in cs.model_sims
        if !(Interfacer.name(sim) == "ConstantIce")
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
function reset_time!(cs::Interfacer.CoupledSimulation, t)
    for sim in cs.model_sims
        Interfacer.reinit!(sim.integrator, sim.integrator.u, t0=t)
    end
end

"""Restarts simulations."""
function restart_sims!(cs::Interfacer.CoupledSimulation)
    @info "Reading checkpoint!"
    t = Dates.datetime2epochms(cs.dates.date[1])
    t0 = Dates.datetime2epochms(cs.dates.date0[1])
    time = Int((t - t0) / 1e3)
    rename_files(cs, 0, time, true)
    Checkpointer.restart!(cs, cs.dirs.checkpoints, time)
    rename_files(cs, 0, time)
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

function update_ocean_values!(cs, ice_T)
    atmos_states = copy(cs.model_sims.atmos_sim.integrator.sol.u)
    atmos_vals = extract_matrix(atmos_states, "atm")
    bound_atmos_vals = atmos_vals[1, :]
    if cs.model_sims.ocean_sim.params.boundary_mapping == "mean"
        atmos_T = mean(bound_atmos_vals)
        update_field!(cs.model_sims.ocean_sim, atmos_T, ice_T)
    else
        update_field!(cs.model_sims.ocean_sim, bound_atmos_vals, ice_T)
    end
    return atmos_vals, bound_atmos_vals
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
    cs.dates.date[1] = TimeManager.current_date(cs, tspan[begin])

    @info("Starting coupling loop")

    # Extract the initial values to be able to stop the simulation if the model goes unstable
    starting_temp_atm = parent(cs.model_sims.atmos_sim.Y_init)
    starting_temp_oce = parent(cs.model_sims.ocean_sim.Y_init)
    starting_temp_ice = parent(cs.model_sims.ice_sim.Y_init)
    upper_limit_temp = maximum([
        maximum(starting_temp_oce),
        maximum(starting_temp_atm),
        maximum(starting_temp_ice),
    ])
    lower_limit_temp = minimum([
        minimum(starting_temp_oce),
        minimum(starting_temp_atm),
        minimum(starting_temp_ice),
    ])

    for t in ((tspan[begin]+Δt_cpl):Δt_cpl:tspan[end])
        time = Int(t - Δt_cpl)
        @info(cs.dates.date[1])

        # Checkpoint to save initial values at this coupling step
        Checkpointer.checkpoint_sims(cs)
        rename_files(cs, 0, time)

        iter = 1
        atmos_vals_list = []
        ocean_vals_list = []
        bound_atmos_vals = nothing
        bound_ocean_vals = nothing
        stopped_at_nan_atm = false
        stopped_at_nan_oce = false

        ice_T = get_field(cs.model_sims.ice_sim, Val(:T_ice))

        while true
            @info("Current iter: $(iter)")

            # Temperature values for the previous iteration.
            pre_bound_atmos_vals = bound_atmos_vals
            pre_bound_ocean_vals = bound_ocean_vals

            if parallel
                FieldExchanger.step_model_sims!(cs.model_sims, t)
                atmos_vals, bound_atmos_vals = update_atmos_values!(cs, ice_T)
                ocean_vals, bound_ocean_vals = update_ocean_values!(cs, ice_T)
            else
                Interfacer.step!(cs.model_sims.ice_sim, t)
                Interfacer.step!(cs.model_sims.ocean_sim, t)
                ocean_vals, bound_ocean_vals = update_atmos_values!(cs, ice_T)

                Interfacer.step!(cs.model_sims.atmos_sim, t)
                atmos_vals, bound_atmos_vals = update_ocean_values!(cs, ice_T)
            end

            stable, stopped_at_nan_atm, stopped_at_nan_oce = is_stable(
                atmos_vals,
                ocean_vals,
                upper_limit_temp,
                lower_limit_temp,
                iter,
            )
            if !stable
                break
            end

            if has_converged(
                bound_atmos_vals,
                pre_bound_atmos_vals,
                bound_ocean_vals,
                pre_bound_ocean_vals,
                iter,
            )
                break
            end

            # Update lists for convergence factor computation
            push!(atmos_vals_list, bound_atmos_vals)
            push!(ocean_vals_list, bound_ocean_vals)

            Checkpointer.checkpoint_sims(cs)
            rename_files(cs, iter, time)

            if iter == iterations
                @info("Stopped at iter $iter due to limit on iterations")
                break
            end
            iter += 1

            # Reset values to beginning of coupling time step
            restart_sims!(cs)
            reset_time!(cs, t - Δt_cpl)
        end
        # Update time and restart integrator at t with the current temperature value
        cs.dates.date[1] = TimeManager.current_date(cs, t)
        if t != tspan[end]
            # TODO: This is a bit sneaky as it is not allowed if t >= sim.integrator.t
            # However, it appears that sim.integrator.t always overshoots and is in fact larger.
            # Is there a better method?
            # It is needed because when running reset_time (reinit) all previous solution values are deleted.
            # In the last iteration, we do not reset, and thus already have values in the solution 
            # for the next time step. So the first iteration solution in this time step have additional 
            # values, if this is not run to empty it.
            reset_time!(cs, t)
        end
        if iterations > 1
            conv_fac_atm, conv_fac_oce = compute_ρ(atmos_vals_list, ocean_vals_list, stopped_at_nan_atm, stopped_at_nan_oce)
            return conv_fac_atm, conv_fac_oce
        end
    end
    return nothing, nothing
end

function compute_ρ(atmos_vals_list, ocean_vals_list, stopped_at_nan_atm, stopped_at_nan_oce)
    if stopped_at_nan_atm || stopped_at_nan_oce
        conv_fac_atm = stopped_at_nan_atm ? Inf : NaN
        conv_fac_oce = stopped_at_nan_oce ? Inf : NaN
        return conv_fac_atm, conv_fac_oce
    end

    conv_fac_atm = []
    conv_fac_oce = []
    bound_error_atm = 0
    bound_error_oce = 0
    end_of_loop = length(atmos_vals_list) - 1
    for i = 1:end_of_loop
        # Compute Schwarz iteration errors
        pre_bound_error_atm = bound_error_atm
        pre_bound_error_oce = bound_error_oce
        bound_error_atm = abs.(atmos_vals_list[i] .- atmos_vals_list[end])
        bound_error_oce = abs.(ocean_vals_list[i] .- ocean_vals_list[end])
        tols_atm =
            100 * eps.(max.(abs.(atmos_vals_list[i]), abs.(atmos_vals_list[end])))
        tols_oce =
            100 * eps.(max.(abs.(ocean_vals_list[i]), abs.(ocean_vals_list[end])))

        # Compute convergence factor
        if i > 1
            indices_atm = findall(
                (pre_bound_error_atm .>= tols_atm) .&
                (bound_error_atm .>= tols_atm),
            )
            indices_oce = findall(
                (pre_bound_error_oce .>= tols_oce) .&
                (pre_bound_error_oce .>= tols_oce),
            )
            conv_fac_atm_value = sqrt(
                sum(bound_error_atm[indices_atm][1:end-1] .^ 2) ./
                sum(pre_bound_error_atm[indices_atm][1:end-1] .^ 2),
            )
            conv_fac_oce_value = sqrt(
                sum(bound_error_oce[indices_oce][1:end-1] .^ 2) ./
                sum(pre_bound_error_oce[indices_oce][1:end-1] .^ 2),
            )
            push!(conv_fac_atm, conv_fac_atm_value)
            push!(conv_fac_oce, conv_fac_oce_value)
        end
    end
    return conv_fac_atm, conv_fac_oce
end

function update_ρ_parallel(conv_fac_atm, conv_fac_oce)
    conv_fac_atm = conv_fac_atm[1:end-1] .* conv_fac_atm[2:end]
    conv_fac_oce = conv_fac_oce[1:end-1] .* conv_fac_oce[2:end]
    conv_fac_atm[1:2:end] .= NaN
    conv_fac_oce[2:2:end] .= NaN
    return conv_fac_atm, conv_fac_oce
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
    params=Dict{Symbol,Int}(),
)
    physical_values = define_realistic_vals()
    merge!(physical_values, params)
    correct_for_a_i!(physical_values)
    compute_C_AO!(physical_values)

    cs = get_coupled_sim(physical_values)
    conv_fac_atm, conv_fac_oce = solve_coupler!(
        cs,
        iterations=iterations,
        parallel=parallel,
    )
    return cs, conv_fac_atm, conv_fac_oce
end;
