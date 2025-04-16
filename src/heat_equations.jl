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

        while true
            @info("Current iter: $(iter)")

            if !parallel
                # Temperature values for the previous iteration.
                pre_bound_atmos_vals = bound_atmos_vals
                pre_bound_ocean_vals = bound_ocean_vals

                # Update ice and ocean models
                Interfacer.step!(cs.model_sims.ice_sim, t) # TODO: Allow for an update
                Interfacer.step!(cs.model_sims.ocean_sim, t)

                # Update ocean value for atmosphere simulation based on boundary mapping
                ice_T = get_field(cs.model_sims.ice_sim, Val(:T_ice))
                ocean_states = copy(cs.model_sims.ocean_sim.integrator.sol.u)
                ocean_vals = extract_matrix(ocean_states, "oce")
                bound_ocean_vals = ocean_vals[end, :]
                ocean_T = mean(bound_ocean_vals)
                if cs.model_sims.atmos_sim.params.boundary_mapping == "mean"
                    update_field!(cs.model_sims.atmos_sim, ocean_T, ice_T)
                else
                    update_field!(cs.model_sims.atmos_sim, bound_ocean_vals, ice_T)
                end

                # Update atmosphere model 
                Interfacer.step!(cs.model_sims.atmos_sim, t)

                # Update atmosphere value for ocean simulation based on boundary mapping
                atmos_states = copy(cs.model_sims.atmos_sim.integrator.sol.u)
                atmos_vals = extract_matrix(atmos_states, "atm")
                bound_atmos_vals = atmos_vals[1, :]
                if cs.model_sims.ocean_sim.params.boundary_mapping == "mean"
                    atmos_T = mean(bound_atmos_vals)
                    update_field!(cs.model_sims.ocean_sim, atmos_T, ice_T)
                else
                    update_field!(cs.model_sims.ocean_sim, bound_atmos_vals, ice_T)
                end

                # If the model goes unstable, abort simulation
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

                # Check for convergence
                if iter > 1
                    converged = has_converged(
                        bound_atmos_vals,
                        pre_bound_atmos_vals,
                        bound_ocean_vals,
                        pre_bound_ocean_vals,
                        iter,
                    )
                    if converged
                        break
                    end
                end

                # Update lists for convergence factor computation
                push!(atmos_vals_list, bound_atmos_vals)
                push!(ocean_vals_list, bound_ocean_vals)

                # Checkpoint to save progress
                Checkpointer.checkpoint_sims(cs)
                rename_files(cs, iter, time)

                # Stop at maximum number of iterations
                if iter == iterations
                    println("Stopped at iter $iter due to limit on iterations")
                    break
                end
                iter += 1

                # Reset the temperature to that of the first checkpoint at this coupling step
                restart_sims!(cs)

                # Reset the integrator time and temperature. 
                reset_time!(cs, t - Δt_cpl)
            else
                # Update all models
                FieldExchanger.step_model_sims!(cs.model_sims, t)

                # Extract atmosphere and ocean boundary temperatures
                atmos_states = copy(cs.model_sims.atmos_sim.integrator.sol.u)
                ocean_states = copy(cs.model_sims.ocean_sim.integrator.sol.u)
                pre_bound_atmos_vals = bound_atmos_vals
                pre_bound_ocean_vals = bound_ocean_vals
                atmos_vals = extract_matrix(atmos_states, "atm")
                ocean_vals = extract_matrix(ocean_states, "oce")
                bound_atmos_vals = atmos_vals[1, :]
                bound_ocean_vals = ocean_vals[end, :]

                # If the model goes unstable, abort simulation
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

                # Check for convergence
                if iter > 1
                    converged = has_converged(
                        bound_atmos_vals,
                        pre_bound_atmos_vals,
                        bound_ocean_vals,
                        pre_bound_ocean_vals,
                        iter,
                    )
                    if converged
                        break
                    end
                end

                # Update lists for convergence factor computation
                push!(atmos_vals_list, bound_atmos_vals)
                push!(ocean_vals_list, bound_ocean_vals)

                # Checkpoint to save progress
                Checkpointer.checkpoint_sims(cs) # I had to remove nothing here
                rename_files(cs, iter, time)

                # Updating the simulations based on the other simulation boundary value.
                # note: only checking the atmosphere mapping here, but the mappings should be the same.
                if cs.model_sims.atmos_sim.params.boundary_mapping == "mean"
                    ocean_T = mean([ocean_state[end] for ocean_state in ocean_states])
                    atmos_T = mean([atmos_state[1] for atmos_state in atmos_states])
                else
                    ocean_T = bound_ocean_vals
                    atmos_T = bound_atmos_vals
                end
                ice_T = get_field(cs.model_sims.ice_sim, Val(:T_ice))
                update_field!(cs.model_sims.ice_sim, Val(:T_ice), ice_T)
                update_field!(cs.model_sims.atmos_sim, ocean_T, ice_T)
                update_field!(cs.model_sims.ocean_sim, atmos_T, ice_T)

                # Stop at maximum number of iterations
                if iter == iterations
                    println("Stopped at iter $iter due to limit on iterations")
                    break
                end
                iter += 1

                # Reset the temperature to that of the first checkpoint at this coupling step
                restart_sims!(cs)

                # Reset the integrator time and temperature. 
                reset_time!(cs, t - Δt_cpl)
            end
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
-`boundary_mapping::String`: Determines if mean or closest in time boundary mapping is used, default: `"mean"`.
-`params::Dict`: Physical parameters and time stepping parameters, default: `false`, default: `Dict{Symbol,Int}()`.
-`analytic_conv_fac::Boolean`: Whether to compute and plot or print the analytical convergence factor, default: `false`.
-`compute_atm_conv_fac::Boolean`: Whether to consider the atmospheric convergence factor, default: `true`.
-`compute_oce_conv_fac::Boolean`: Whether to consider the oceanic convergence factor, default: `true`.
-`plot_unstable_range::Boolean`: Whether to plot the unstable region (cfl condition), default: `false`.
-`a_is::Array`: If not empty, the convergence factor is plotted as a function of `a_is`, default: `[]`.
-`var_name::String`: If not nothing, the convergence factor is plotted wrt the variable.
    If there are also several `a_is`, it is plotted for all of the `a_i` in the same plot, default: `nothing`.
-`xscale::Symbol`: Plotting argument for x-scale, default: `:identity`.
-`yscale::Symbol`: Plotting argument for y-scale, default: `:identity`.
-`legend::Symbol`: Legend position, default: `:right`.
-`xticks::Symbol or Array`: Defines the tick marks on the x-axis, default: `nothing`.
-`yticks::Symbol or Array`: Defines the tick marks on the y-axis, default: `:auto`.
-`text_scaling::Tuple`: If the model goes unstable, there is text in the plot. This argument 
    changes the text position, default: `(1, 5)`.

"""
function coupled_heat_equations(;
    iterations=1,
    parallel=false,
    boundary_mapping="mean",
    params=Dict{Symbol,Int}(),
    analytic_conv_fac=false,
    compute_atm_conv_fac=true,
    compute_oce_conv_fac=true,
    a_is=[],
    var_name=nothing,
    xscale=:identity,
    yscale=:identity,
    legend=:right,
    xticks=nothing,
    yticks=:auto,
    text_scaling=(1, 5),
)
    physical_values = define_realistic_vals()
    merge!(physical_values, params)
    correct_for_a_i!(physical_values)
    compute_C_AO!(physical_values)
    physical_values[:boundary_mapping] = boundary_mapping

    if !(
        !isempty(a_is) ||
        !isnothing(var_name)
    )
        # Run coupled simulation
        cs = get_coupled_sim(NamedTuple(physical_values))
        conv_fac_atm, conv_fac_oce = solve_coupler!(
            cs,
            iterations=iterations,
            parallel=parallel,
        )
        return cs, conv_fac_atm, conv_fac_oce

    elseif !isnothing(var_name)
        # Plot convergence factor with respect to some parameter, and different a_i
        variable_dict = get_var_dict()
        color_dict, linestyle_dict = get_color_dict()
        var = variable_dict[Symbol(var_name)][1]
        conv_facs_atm, conv_facs_oce, param_analytic, conv_facs_analytic =
            isempty(a_is) ?
            get_conv_facs_one_variable(
                physical_values,
                var,
                var_name,
                iterations=iterations,
                analytic=analytic_conv_fac,
                log_scale=(xscale == :log10),
            ) :
            get_conv_facs_one_variable(
                physical_values,
                var,
                var_name,
                iterations=iterations,
                a_i_variable=a_is,
                analytic=analytic_conv_fac,
                log_scale=(xscale == :log10),
            )
        var, conv_facs_oce, conv_facs_atm, param_analytic, conv_facs_analytic =
            handle_variable(
                var,
                var_name,
                conv_facs_oce,
                conv_facs_atm,
                physical_values;
                dims=2,
                param_analytic=param_analytic,
                conv_facs_analytic=conv_facs_analytic,
            )

        xticks = !isnothing(xticks) ? xticks : var
        if isempty(a_is)
            plot_wrt_a_i_and_one_param(
                conv_facs_oce,
                conv_facs_atm,
                [physical_values[:a_i]],
                var,
                variable_dict[Symbol(var_name)][2],
                conv_facs_analytic=conv_facs_analytic,
                param_analytic=param_analytic,
                xticks=xticks,
                yticks=yticks,
                xscale=xscale,
                yscale=yscale,
                colors=[color_dict[round(physical_values[:a_i], digits=1)]],
                linestyles=[linestyle_dict[round(physical_values[:a_i], digits=1)]],
                text_scaling=text_scaling,
                legend=legend,
                compute_atm_conv_fac=compute_atm_conv_fac,
                compute_oce_conv_fac=compute_oce_conv_fac,
            )
        else
            plot_wrt_a_i_and_one_param(
                conv_facs_oce,
                conv_facs_atm,
                a_is,
                var,
                variable_dict[Symbol(var_name)][2],
                conv_facs_analytic=conv_facs_analytic,
                param_analytic=param_analytic,
                xticks=xticks,
                yticks=yticks,
                xscale=xscale,
                yscale=yscale,
                colors=[color_dict[round(a_i, digits=1)] for a_i in a_is],
                linestyles=[linestyle_dict[round(a_i, digits=1)] for a_i in a_is],
                text_scaling=text_scaling,
                legend=legend,
                compute_atm_conv_fac=compute_atm_conv_fac,
                compute_oce_conv_fac=compute_oce_conv_fac,
            )
        end
    end
end;






