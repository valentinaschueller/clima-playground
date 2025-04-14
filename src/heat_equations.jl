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
-`print_conv::Boolean`: Whether to print the convergence factor or not, default: `false`.
-`plot_conv::Boolean`: Whether to plot the convergence factor or not, default: `false`.
-`return_conv::Boolean`: Whether to return the convergence factor or not, default: `false`.
-`analytic_conv_fac_value::Float64`: The analytical convergence factor. If sent in,
    it is plotted and printed with the numerical convergence factor, default: `nothing`.
-`combine_ρ_parallel::Boolean`: Whether to combine two on eachother following convergence factor values
    for the parallel Schwarz iteration, default: `false`.
-`compute_atm_conv_fac::Boolean`: Whether to consider the atmospheric convergence factor or only the oceanic, default: `true`.
-`compute_oce_conv_fac::Boolean`: Whether to consider the oceanic convergence factor, default: `true`.
-`legend::Symbol`: Legend placement for plotting, default: `:right`.

"""
function solve_coupler!(
    cs::Interfacer.CoupledSimulation;
    iterations=1,
    parallel=false,
    print_conv=false,
    plot_conv=false,
    return_conv=false,
    analytic_conv_fac_value=nothing,
    combine_ρ_parallel=false,
    compute_atm_conv_fac=true,
    compute_oce_conv_fac=true,
    legend=:right,
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

        # If we want to compute the convergence factor
        if print_conv || plot_conv || return_conv
            # If we had instabilities in one of the models, we cannot compute it
            if stopped_at_nan_atm && stopped_at_nan_oce
                return Inf, Inf
            elseif stopped_at_nan_oce
                return NaN, Inf
            elseif stopped_at_nan_atm
                return Inf, NaN
            else
                end_of_loop = length(atmos_vals_list) - 1
            end

            conv_fac_atm = []
            conv_fac_oce = []
            bound_error_atm = 0
            bound_error_oce = 0
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
            # If we run the parallel iteration and want to plot with the analytical convergence factor,
            # Two on eachother following factors should be combined, this can also be chosen via
            # the combine_ρ_parallel argument.
            if parallel && (!isnothing(analytic_conv_fac_value) || combine_ρ_parallel)
                conv_fac_atm = conv_fac_atm[1:end-1] .* conv_fac_atm[2:end]
                conv_fac_oce = conv_fac_oce[1:end-1] .* conv_fac_oce[2:end]
                conv_fac_atm[1:2:end] .= NaN
                conv_fac_oce[2:2:end] .= NaN
                if !isnothing(analytic_conv_fac_value)
                    ylabel = L"$\rho_{k+1}\times\rho_k,$ $\hat{\rho}$"
                else
                    ylabel = L"$\rho_{k+1}\times\rho_k$"
                end
            elseif !isnothing(analytic_conv_fac_value)
                ylabel = L"$\hat{\rho}$, $\rho_k$"
            else
                ylabel = L"$\rho_k$"
            end
            if print_conv
                if (!isnothing(analytic_conv_fac_value) || combine_ρ_parallel) && parallel
                    if compute_atm_conv_fac && compute_oce_conv_fac
                        println("Combined convergence factor atmosphere: $conv_fac_atm")
                        println("Combined convergence factor ocean: $conv_fac_oce")
                    elseif compute_oce_conv_fac
                        println("Combined convergence factor ocean: $conv_fac_oce")
                    elseif compute_atm_conv_fac
                        println("Combined convergence factor atmosphere: $conv_fac_atm")
                    end
                else
                    if compute_atm_conv_fac && compute_oce_conv_fac
                        println("Convergence factor ocean: $conv_fac_oce")
                        println("Convergence factor atmosphere: $conv_fac_atm")
                    elseif compute_oce_conv_fac
                        println("Convergence factor ocean: $conv_fac_oce")
                    elseif compute_atm_conv_fac
                        println("Convergence factor atmosphere: $conv_fac_atm")
                    end
                end
                if !isnothing(analytic_conv_fac_value) && !parallel
                    println("Analytical convergence factor: $analytic_conv_fac_value")
                end
            end
            if plot_conv
                gr()
                plot()
                color_dict, _ = get_color_dict()
                color = color_dict[round(cs.model_sims.atmos_sim.params.a_i, digits=1)]
                k_atm = 2:length(conv_fac_atm)+1
                k_oce = 2:length(conv_fac_oce)+1
                if compute_atm_conv_fac && compute_oce_conv_fac
                    scatter!(
                        k_atm,
                        conv_fac_atm,
                        label="atm",
                        legend=legend,
                        color=color,
                        markershape=:x,
                        markersize=5,
                        xlabel="k",
                        ylabel=ylabel,
                        ylim=(
                            0,
                            maximum([
                                maximum(conv_fac_atm[.!isnan.(conv_fac_atm)]),
                                maximum(conv_fac_oce[.!isnan.(conv_fac_oce)]),
                            ]) * 1.2,
                        ),
                    )
                    scatter!(
                        k_oce,
                        conv_fac_oce,
                        label="oce",
                        legend=legend,
                        color=color,
                        markershape=:circle,
                        markersize=5,
                    )
                elseif compute_oce_conv_fac
                    scatter!(
                        k_oce,
                        conv_fac_oce,
                        label="oce",
                        legend=legend,
                        color=color,
                        markershape=:circle,
                        markersize=5,
                        xlabel="k",
                        ylabel=ylabel,
                        ylim=(0, maximum(conv_fac_oce[.!isnan.(conv_fac_oce)]) * 1.2),
                    )
                elseif compute_atm_conv_fac
                    scatter!(
                        k_atm,
                        conv_fac_atm,
                        label="atm",
                        legend=legend,
                        color=color,
                        markershape=:circle,
                        markersize=5,
                        xlabel="k",
                        ylabel=ylabel,
                        ylim=(0, maximum(conv_fac_atm[.!isnan.(conv_fac_atm)]) * 1.2),
                    )
                end
                if !isnothing(analytic_conv_fac_value)
                    scatter!(
                        k_atm,
                        ones(length(k_atm)) * analytic_conv_fac_value,
                        label="analytic",
                        color=color,
                        markershape=:hline,
                        ylim=(0, analytic_conv_fac_value * 1.2),
                    )
                end
                display(current())

            end
            if return_conv
                if compute_atm_conv_fac && compute_oce_conv_fac
                    return conv_fac_atm, conv_fac_oce
                elseif compute_oce_conv_fac
                    return nothing, conv_fac_oce
                elseif compute_atm_conv_fac
                    return conv_fac_atm, nothing
                end
            end
        end
    end
end

"""
Setup for running the coupled simulation and running it.

**Optional Keyword Arguments:**

-`iterations::Int`: Number of iterations before the Schwarz iteration is terminated, default: `1`.
-`parallel::Boolean`: Whether to run the parallel or alternating Schwarz iteration, default: `false`.
-`boundary_mapping::String`: Determines if mean or closest in time boundary mapping is used, default: `"mean"`.
-`params::Dict`: Physical parameters and time stepping parameters, default: `false`, default: `Dict{Symbol,Int}()`.
-`print_conv_facs_iter::Boolean`: Whether to print the convergence factor or not, default: `false`.
-`plot_conv_facs_iter::Boolean`: Whether to plot the convergence factor or not, default: `false`.
-`analytic_conv_fac::Boolean`: Whether to compute and plot or print the analytical convergence factor, default: `false`.
-`compute_atm_conv_fac::Boolean`: Whether to consider the atmospheric convergence factor, default: `true`.
-`compute_oce_conv_fac::Boolean`: Whether to consider the oceanic convergence factor, default: `true`.
-`combine_ρ_parallel::Boolean`: Whether to combine two on eachother following convergence factor values
    for the parallel Schwarz iteration, default: `false`.
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
    print_conv_facs_iter=false,
    plot_conv_facs_iter=false,
    analytic_conv_fac=false,
    compute_atm_conv_fac=true,
    compute_oce_conv_fac=true,
    combine_ρ_parallel=false,
    plot_unstable_range=false,
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
    if !isempty(params)
        merge!(physical_values, params)
        if haskey(params, :a_i)
            physical_values = update_physical_values(params[:a_i], physical_values)
        end
    end
    physical_values[:boundary_mapping] = boundary_mapping

    if !(
        plot_conv_facs_iter ||
        print_conv_facs_iter ||
        plot_unstable_range ||
        !isempty(a_is) ||
        !isnothing(var_name)
    )
        # Run coupled simulation
        cs = get_coupled_sim(
            physical_values,
            boundary_mapping=physical_values[:boundary_mapping],
        )
        solve_coupler!(
            cs,
            iterations=iterations,
            parallel=parallel,
            compute_atm_conv_fac=compute_atm_conv_fac,
            compute_oce_conv_fac=compute_oce_conv_fac,
        )
        if analytic_conv_fac
            println(compute_ρ_analytical(physical_values))
        end

    elseif plot_conv_facs_iter || print_conv_facs_iter
        # Compute and plot or print convergence factor with respect to iteration
        cs = get_coupled_sim(
            physical_values,
            boundary_mapping=physical_values[:boundary_mapping],
        )
        if analytic_conv_fac
            analytic_conv_fac_value = compute_ρ_analytical(physical_values)
        else
            analytic_conv_fac_value = nothing
        end
        solve_coupler!(
            cs,
            parallel=parallel,
            iterations=iterations,
            plot_conv=plot_conv_facs_iter,
            print_conv=print_conv_facs_iter,
            analytic_conv_fac_value=analytic_conv_fac_value,
            combine_ρ_parallel=combine_ρ_parallel,
            compute_atm_conv_fac=compute_atm_conv_fac,
            compute_oce_conv_fac=compute_oce_conv_fac,
            legend=legend,
        )

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

    elseif !isempty(a_is) && !plot_unstable_range
        # Plot convergence factor wrt a_i
        conv_facs_atm, conv_facs_oce, param_analytic, conv_facs_analytic =
            get_conv_facs_one_variable(
                physical_values,
                a_is,
                "a_i",
                iterations=iterations,
                analytic=analytic_conv_fac,
                log_scale=(xscale == :log10),
            )
        xticks = !isnothing(xticks) ? xticks : a_is
        plot_wrt_a_i_and_one_param(
            conv_facs_oce,
            conv_facs_atm,
            [physical_values[:a_i]],
            a_is,
            L"$a^I$",
            conv_facs_analytic=conv_facs_analytic,
            param_analytic=param_analytic,
            xticks=xticks,
            yticks=yticks,
            xscale=xscale,
            yscale=yscale,
            colors=[:black],
            linestyles=[:solid],
            text_scaling=text_scaling,
            legend=legend,
            compute_atm_conv_fac=compute_atm_conv_fac,
            compute_oce_conv_fac=compute_oce_conv_fac,
        )

    elseif plot_unstable_range
        # Plot CFL condition range for atmosphere or ocean
        variable_dict = get_var_dict()
        color_dict, linestyle_dict = get_color_dict()
        physical_values[:Δt_min] = 100

        Δz = 10 .^ LinRange(log10(0.001), log10(1), 50)
        Δt = 10 .^ LinRange(log10(1), log10(100), 50)
        n_zs_atm =
            compute_atm_conv_fac ?
            Int.(
                round.(
                    (physical_values[:h_atm] - physical_values[:z_0numA]) ./ reverse(Δz)
                )
            ) : nothing
        n_ts_atm =
            compute_atm_conv_fac ? Int.(round.(physical_values[:Δt_min] ./ reverse(Δt))) :
            nothing
        n_zs_oce =
            compute_oce_conv_fac ?
            Int.(
                round.(
                    (physical_values[:h_oce] - physical_values[:z_0numO]) ./ reverse(Δz)
                )
            ) : nothing
        n_ts_oce =
            compute_oce_conv_fac ? Int.(round.(physical_values[:Δt_min] ./ reverse(Δt))) :
            nothing

        xscale = :log10
        yscale = :log10
        xticks = [1, 10, 100]
        yticks = [0.001, 0.01, 0.1]

        a_is = !isempty(a_is) ? a_is : [physical_values[:a_i]]

        plot()
        for a_i in a_is
            physical_values[:a_i] = a_i
            if compute_atm_conv_fac
                unstable_matrix_atm =
                    stability_check(physical_values, n_zs_atm, n_ts_atm, "n_atm", "n_t_atm")
                Δz, _, unstable_matrix_atm, _, _ = handle_variable(
                    n_zs_atm,
                    "n_atm",
                    nothing,
                    unstable_matrix_atm,
                    physical_values;
                    dims=1,
                )
                Δt, _, unstable_matrix_atm, _, _ = handle_variable(
                    n_ts_atm,
                    "n_t_atm",
                    nothing,
                    unstable_matrix_atm,
                    physical_values;
                    dims=2,
                )

                plot_Δz_Δt(
                    unstable_matrix_atm,
                    Δz,
                    Δt,
                    L"$\Delta z^A$",
                    L"$\Delta t^A$",
                    xscale=xscale,
                    yscale=yscale,
                    xticks=xticks,
                    yticks=yticks,
                    color=color_dict[round(a_i, digits=1)],
                    a_i=a_i,
                    legend=legend,
                )
            elseif compute_oce_conv_fac
                unstable_matrix_oce =
                    stability_check(physical_values, n_zs_oce, n_ts_oce, "n_oce", "n_t_oce")
                Δz, unstable_matrix_oce, _, _, _ = handle_variable(
                    n_zs_oce,
                    "n_oce",
                    unstable_matrix_oce,
                    nothing,
                    physical_values;
                    dims=1,
                )
                Δt, unstable_matrix_oce, _, _, _ = handle_variable(
                    n_ts_oce,
                    "n_t_oce",
                    unstable_matrix_oce,
                    nothing,
                    physical_values;
                    dims=2,
                )
                plot_Δz_Δt(
                    unstable_matrix_oce,
                    Δz,
                    Δt,
                    L"$\Delta z^O$",
                    L"$\Delta t^O$",
                    xscale=xscale,
                    yscale=yscale,
                    xticks=xticks,
                    yticks=yticks,
                    color=color_dict[round(a_i, digits=1)],
                    a_i=a_i,
                    legend=legend,
                )
            end
        end
        display(current())
    end
end;






