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
function reinit!(cs::Interfacer.CoupledSimulation, t)
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

function set_time!(cs::Interfacer.CoupledSimulation, t)
    cs.dates.date[1] = TimeManager.current_date(cs, t)
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

    ρ_atm = nothing
    ρ_oce = nothing

    for t in ((tspan[begin]+Δt_cpl):Δt_cpl:tspan[end])
        time = Int(t - Δt_cpl)
        set_time!(cs, time)
        @info("Current time: $time")

        # Sets u0 = u(time), t0 = time, and empties u
        reinit!(cs, time)

        # Checkpoint to save initial values at this coupling step
        Checkpointer.checkpoint_sims(cs)
        rename_files(cs, 0, time)

        iter = 1
        atmos_vals_list = []
        ocean_vals_list = []
        bound_atmos_vals = nothing
        bound_ocean_vals = nothing

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

            if (
                !is_stable(atmos_vals, upper_limit_temp, lower_limit_temp) ||
                !is_stable(ocean_vals, upper_limit_temp, lower_limit_temp)
            )
                break
            end

            push!(atmos_vals_list, bound_atmos_vals)
            push!(ocean_vals_list, bound_ocean_vals)

            if has_converged(
                bound_atmos_vals,
                pre_bound_atmos_vals,
                bound_ocean_vals,
                pre_bound_ocean_vals,
                iter,
            )
                break
            end

            Checkpointer.checkpoint_sims(cs)
            rename_files(cs, iter, time)

            if iter == iterations
                @info("Stopped at iter $iter due to limit on iterations")
                break
            end
            iter += 1

            # Reset values to beginning of coupling time step
            restart_sims!(cs)
            reinit!(cs, time)
        end
        if iterations > 1
            ρ_atm, ρ_oce = compute_ρ_numerical(atmos_vals_list, ocean_vals_list)
        end
    end
    set_time!(cs, tspan[end])
    return ρ_atm, ρ_oce
end

function run_simulation(
    physical_values;
    iterations=10,
    parallel=false,
)
    cs = get_coupled_sim(physical_values)
    ρ_atm, ρ_oce = solve_coupler!(
        cs,
        iterations=iterations,
        parallel=parallel,
    )
    return cs, ρ_atm, ρ_oce
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

    cs, ρ_atm, ρ_oce = run_simulation(physical_values, iterations=iterations, parallel=parallel)
    return cs, ρ_atm, ρ_oce
end;
