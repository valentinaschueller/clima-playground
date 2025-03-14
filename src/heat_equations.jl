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
    Checkpointer,
    FieldExchanger,
    FluxCalculator,
    Interfacer,
    TimeManager,
    Utilities

function rename_file(cs::Interfacer.CoupledSimulation, iter, time, reverse=false)
    """When a file is saved, its always called the same thing.
    Had to rename it for each iteration to not overwrite"""

    for sim in cs.model_sims
        if !(Interfacer.name(sim) == "ConstantIce")
            original_file = joinpath(cs.dirs.artifacts, "checkpoints", "checkpoint_" * Interfacer.name(sim) * "_$time.hdf5")
            new_file = joinpath(cs.dirs.artifacts, "checkpoints", "checkpoint_" * Interfacer.name(sim) * "_$iter" * "_$time.hdf5")
            if !reverse
                mv(original_file, new_file, force=true)
            else
                mv(new_file, original_file, force=true)
            end
        end
    end
end

function reset_time!(cs::Interfacer.CoupledSimulation, t)
    """resets integrator time"""
    for sim in cs.model_sims
        Interfacer.reinit!(sim.integrator, sim.integrator.u, t0=t)
    end
end

function restart_sims!(cs::Interfacer.CoupledSimulation)
    @info "Reading checkpoint!"
    t = Dates.datetime2epochms(cs.dates.date[1])
    t0 = Dates.datetime2epochms(cs.dates.date0[1])
    time = Int((t - t0) / 1e3)
    for sim in cs.model_sims
        if Checkpointer.get_model_prog_state(sim) !== nothing
            rename_file(cs, 0, time, true)
            Checkpointer.restart_model_state!(sim, cs.comms_ctx, time, input_dir=cs.dirs.artifacts)
            rename_file(cs, 0, time)
        end
    end
end

function solve_coupler!(cs::Interfacer.CoupledSimulation; iterate=10, parallel=false, print_conv=false, plot_conv=false, return_conv=false, analytic_conv_fac_value=nothing, combine=false, atm=true, oce=true, legend=:topright)
    (; Δt_cpl, tspan) = cs
    """
    Runs the coupled CoupledSimulation

    cs: Interfacer.CoupledSimulation, A coupled simulation with atmosphere, ice and ocean

    iterate: Int, Number of iterations before the Schwarz iteration is terminated
    parallel: boolean, Whether to run the parallel or alternating Schwarz iteration
    print_conv: boolean, Whether to print the convergence factor or not
    plot_conv: boolean, Whether to plot the convergence factor or not
    return_conv: boolean, Whether to return the convergence factor or not
    analytic_conv_fac_value: Float64, the analytical convergence factor. If sent in,
        it is plotted and printed with the numerical convergence factor
    combine: boolean, Whether to combine two on eachother following convergence factor values
        for the parallel Schwarz iteration
    atm:boolean, Whether to consider the atmospheric convergence factor or only the oceanic
    oce:boolean, Whether to consider the oceanic convergence factor
    """

    cs.dates.date[1] = TimeManager.current_date(cs, tspan[begin])

    @info("Starting coupling loop")

    # Extract the initial values to be able to stop the simulation if the model goes unstable
    starting_temp_oce = cs.model_sims.atmos_sim.params.T_atm_ini
    starting_temp_atm = cs.model_sims.ocean_sim.params.T_oce_ini
    starting_temp_ice = cs.model_sims.atmos_sim.params.T_ice_ini
    upper_limit_temp = maximum([starting_temp_oce, starting_temp_atm, starting_temp_ice])
    lower_limit_temp = minimum([starting_temp_oce, starting_temp_atm, starting_temp_ice])

    for t in ((tspan[begin]+Δt_cpl):Δt_cpl:tspan[end])
        time = Int(t - Δt_cpl)
        @info(cs.dates.date[1])

        # Checkpoint to save initial values at this coupling step
        Checkpointer.checkpoint_sims(cs)
        rename_file(cs, 0, time)

        iter = 1
        atmos_vals_list = []
        ocean_vals_list = []
        bound_atmos_vals = nothing
        bound_ocean_vals = nothing
        stopped_at_nan_atm = false
        stopped_at_nan_oce = false

        # Should handle when it doesnt converge here (for example when delta z too small.)
        while true
            @info("Current iter: $(iter)")

            # Alternating or Parallel Schwarz
            if !parallel
                # Temperature values for the previous iteration.
                pre_bound_atmos_vals = bound_atmos_vals
                pre_bound_ocean_vals = bound_ocean_vals

                # Update ica and ocean models
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
                stable, stopped_at_nan_atm, stopped_at_nan_oce = is_stable(atmos_vals, ocean_vals, upper_limit_temp, lower_limit_temp, iter)
                if !stable
                    break
                end

                # Check for convergence
                if iter > 1
                    converged = has_converged(bound_atmos_vals, pre_bound_atmos_vals, bound_ocean_vals, pre_bound_ocean_vals, iter)
                    if converged
                        break
                    end
                end

                # Update lists for convergence factor computation
                push!(atmos_vals_list, bound_atmos_vals)
                push!(ocean_vals_list, bound_ocean_vals)

                # Checkpoint to save progress
                Checkpointer.checkpoint_sims(cs)
                rename_file(cs, iter, time)

                # Stop at iterate iterations
                if iter == iterate
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

                # Extract atmosphere and ocean boundary
                atmos_states = copy(cs.model_sims.atmos_sim.integrator.sol.u)
                ocean_states = copy(cs.model_sims.ocean_sim.integrator.sol.u)
                pre_bound_atmos_vals = bound_atmos_vals
                pre_bound_ocean_vals = bound_ocean_vals
                atmos_vals = extract_matrix(atmos_states, "atm")
                ocean_vals = extract_matrix(ocean_states, "oce")
                bound_atmos_vals = atmos_vals[1, :]
                bound_ocean_vals = ocean_vals[end, :]

                # If the model goes unstable, abort simulation
                stable, stopped_at_nan_atm, stopped_at_nan_oce = is_stable(atmos_vals, ocean_vals, upper_limit_temp, lower_limit_temp, iter)
                if !stable
                    break
                end

                # Check for convergence
                if iter > 1
                    converged = has_converged(bound_atmos_vals, pre_bound_atmos_vals, bound_ocean_vals, pre_bound_ocean_vals, iter)
                    if converged
                        break
                    end
                end

                # Update lists for convergence factor computation
                push!(atmos_vals_list, bound_atmos_vals)
                push!(ocean_vals_list, bound_ocean_vals)

                # Checkpoint to save progress
                Checkpointer.checkpoint_sims(cs) # I had to remove nothing here
                rename_file(cs, iter, time)

                # Updating the simulations based on the other simulation boundary value.
                # OBS: Only checking the atmosphere mapping here, but the mappings should be the same.
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

                # Stop at iterate iterations
                if iter == iterate
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
                tols_atm = 100 * eps.(max.(abs.(atmos_vals_list[i]), abs.(atmos_vals_list[end])))
                tols_oce = 100 * eps.(max.(abs.(ocean_vals_list[i]), abs.(ocean_vals_list[end])))

                # Compute convergence factor
                if i > 1
                    indices_atm = findall((pre_bound_error_atm .>= tols_atm) .& (bound_error_atm .>= tols_atm))
                    indices_oce = findall((pre_bound_error_oce .>= tols_oce) .& (pre_bound_error_oce .>= tols_oce))
                    conv_fac_atm_value = sum(bound_error_atm[indices_atm]) ./ sum(pre_bound_error_atm[indices_atm])
                    conv_fac_oce_value = sum(bound_error_oce[indices_oce]) ./ sum(pre_bound_error_oce[indices_oce])
                    push!(conv_fac_atm, conv_fac_atm_value)
                    push!(conv_fac_oce, conv_fac_oce_value)
                end
            end
            # If we run the parallel iteration and want to plot with the analytical convergence factor,
            # Two on eachother following factors should be combined, this can also be chosen via
            # the combine argument.
            if parallel && (!isnothing(analytic_conv_fac_value) || combine)
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
                if (!isnothing(analytic_conv_fac_value) || combine) && parallel
                    if atm && oce
                        println("Combined convergence factor atmosphere: $conv_fac_atm")
                        println("Combined convergence factor ocean: $conv_fac_oce")
                    elseif oce
                        println("Combined convergence factor ocean: $conv_fac_oce")
                    elseif atm
                        println("Combined convergence factor atmosphere: $conv_fac_atm")
                    end
                else
                    if atm && oce
                        println("Convergence factor ocean: $conv_fac_oce")
                        println("Convergence factor atmosphere: $conv_fac_atm")
                    elseif oce
                        println("Convergence factor ocean: $conv_fac_oce")
                    elseif atm
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
                if atm && oce
                    scatter!(k_atm, conv_fac_atm, label="atm", legend=legend, color=color, markershape=:x, markersize=5, xlabel="k", ylabel=ylabel, ylim=(0, maximum([maximum(conv_fac_atm[.!isnan.(conv_fac_atm)]), maximum(conv_fac_oce[.!isnan.(conv_fac_oce)])]) * 1.2))
                    scatter!(k_oce, conv_fac_oce, label="oce", legend=legend, color=color, markershape=:circle, markersize=5)
                elseif oce
                    scatter!(k_oce, conv_fac_oce, label="oce", legend=legend, color=color, markershape=:circle, markersize=5, xlabel="k", ylabel=ylabel, ylim=(0, maximum(conv_fac_oce[.!isnan.(conv_fac_oce)]) * 1.2))
                elseif atm
                    scatter!(k_atm, conv_fac_atm, label="atm", legend=legend, color=color, markershape=:circle, markersize=5, xlabel="k", ylabel=ylabel, ylim=(0, maximum(conv_fac_atm[.!isnan.(conv_fac_atm)]) * 1.2))
                end
                if !isnothing(analytic_conv_fac_value)
                    scatter!(k_atm, ones(length(k_atm)) * analytic_conv_fac_value, label="analytic", color=color, markershape=:hline, ylim=(0, analytic_conv_fac_value * 1.2))
                end
                display(current())

            end
            if return_conv
                if atm && oce
                    return conv_fac_atm, conv_fac_oce
                elseif oce
                    return nothing, conv_fac_oce
                elseif atm
                    return conv_fac_atm, nothing
                end
            end
        end
    end
end

# Update such that i can plot both atmosphere and ocean separately, via an oce boolean
# argument. Fix the unstable plot, should be easier now when instability is defined for ocean and atmosphere separately.
# Use is_stable in the unstable method. Then it should be done. Some additional comments and
# Test with all combinations. Perhaps put things in new folders as well.
function coupled_heat_equations(; iterate=10, parallel=false, boundary_mapping="mean", values=Dict{Symbol,Int}(), print_conv_facs_iter=false, plot_conv_facs_iter=false, analytic_conv_fac=false, atm=true, oce=true, combine=false, plot_unstable_range=false, a_is=[], var_name=nothing, xscale=:identity, legend=:topright, log_conv_fac=false, xticks=nothing, text_scaling=(1, 5))
    """
    Setup for running the coupled simulation and running it

    iterate : Int, Number of iterations before the Schwarz iteration is terminated
    parallel : boolean, Whether to run the parallel or alternating Schwarz iteration
    boundary_mapping : String, Determines if mean or closest in time boundary mapping is used
    values : Dict, Physical parameters and time stepping parameters
    print_conv_facs_iter : boolean, Whether to print the convergence factor or not
    plot_conv_facs_iter : boolean, Whether to plot the convergence factor or not
    analytic_conv_fac : boolean, Whether to compute and plot or print the analytical convergence factor.
    atm : boolean, Whether to consider the atmospheric convergence factor
    oce : boolean, Whether to consider the oceanic convergence factor
    combine : boolean, Whether to combine two on eachother following convergence factor values
        for the parallel Schwarz iteration
    plot_unstable_range : boolean, Whether to plot the unstable range (cfl condition)
    a_is : list, If not empty, the convergence factor is plotted as a function of a_is
    var_name : String, If not nothing, the convergence factor is plotted wrt the variable.
        If there are also several a_is, it is plotted for all of the a_i in the same plot
    xscale : Symbol, Plotting argument for x-scale
    legend : Symbol, Legend position
    log_conv_fac: boolean, if the logarithm of the convergence factor should be taken
    text_scaling: tuple, if the model goes unstable, there is text in the plot. This argument 
        changes text position
    """
    # Get physical values
    physical_values = define_realistic_vals()
    if !isempty(values)
        merge!(physical_values, values)
        if haskey(values, :a_i)
            physical_values = update_physical_values(values[:a_i], physical_values)
        end
    end
    physical_values[:boundary_mapping] = boundary_mapping

    if !(plot_conv_facs_iter || print_conv_facs_iter || plot_unstable_range || !isempty(a_is) || !isnothing(var_name))
        # Run coupled simulation
        cs = get_coupled_sim(physical_values, boundary_mapping=physical_values[:boundary_mapping])
        solve_coupler!(cs, iterate=iterate, parallel=parallel, atm=atm, oce=oce)
        if analytic_conv_fac
            println(compute_actual_rho(physical_values))
        end

    elseif plot_conv_facs_iter || print_conv_facs_iter
        # Compute and plot or print convergence factor with respect to iteration
        cs = get_coupled_sim(physical_values, boundary_mapping=physical_values[:boundary_mapping])
        if analytic_conv_fac
            analytic_conv_fac_value = compute_actual_rho(physical_values)
        else
            analytic_conv_fac_value = nothing
        end
        solve_coupler!(cs, parallel=parallel, iterate=iterate, plot_conv=plot_conv_facs_iter, print_conv=print_conv_facs_iter, analytic_conv_fac_value=analytic_conv_fac_value, combine=combine, atm=atm, oce=oce, legend=legend)

    elseif !isnothing(var_name)
        # Plot convergence factor with respect to some parameter, and different a_i
        variable_dict = get_var_dict()
        color_dict, linestyle_dict = get_color_dict()
        var = variable_dict[Symbol(var_name)][1]
        conv_facs_atm, conv_facs_oce, param_analytic, conv_facs_analytic = isempty(a_is) ? get_conv_facs_one_variable(physical_values, var, var_name, analytic=analytic_conv_fac, log_scale=(xscale == :log10)) : get_conv_facs_one_variable(physical_values, var, var_name, a_i_variable=a_is, analytic=analytic_conv_fac, log_scale=(xscale == :log10))
        var, conv_facs_oce, conv_facs_atm, param_analytic, conv_facs_analytic = handle_variable(
            var, var_name, conv_facs_oce, conv_facs_atm, physical_values; dims=2,
            param_analytic=param_analytic, conv_facs_analytic=conv_facs_analytic
        )

        xticks = !isnothing(xticks) ? xticks : var
        if isempty(a_is)
            plot_wrt_a_i_and_one_param(conv_facs_oce, conv_facs_atm, [physical_values[:a_i]], var, variable_dict[Symbol(var_name)][2],
                conv_facs_analytic=conv_facs_analytic, param_analytic=param_analytic, xticks=xticks, log_conv_fac=log_conv_fac, xscale=xscale,
                colors=[color_dict[round(physical_values[:a_i], digits=1)]], linestyles=[linestyle_dict[round(physical_values[:a_i], digits=1)]],
                text_scaling=text_scaling, legend=legend, atm=atm, oce=oce
            )
        else
            plot_wrt_a_i_and_one_param(conv_facs_oce, conv_facs_atm, a_is, var, variable_dict[Symbol(var_name)][2],
                conv_facs_analytic=conv_facs_analytic, param_analytic=param_analytic, xticks=xticks, log_conv_fac=log_conv_fac, xscale=xscale,
                colors=[color_dict[round(a_i, digits=1)] for a_i in a_is], linestyles=[linestyle_dict[round(a_i, digits=1)] for a_i in a_is],
                text_scaling=text_scaling, legend=legend, atm=atm, oce=oce
            )
        end

    elseif !isempty(a_is) && !plot_unstable_range
        # Plot convergence factor wrt a_i
        conv_facs_atm, conv_facs_oce, param_analytic, conv_facs_analytic = get_conv_facs_one_variable(physical_values, a_is, "a_i", analytic=analytic_conv_fac, log_scale=(xscale == :log10))
        xticks = !isnothing(xticks) ? xticks : a_is
        plot_wrt_a_i_and_one_param(conv_facs_oce, conv_facs_atm, [physical_values[:a_i]], a_is, L"$a^I$",
            conv_facs_analytic=conv_facs_analytic, param_analytic=param_analytic, xticks=xticks, log_conv_fac=log_conv_fac, xscale=xscale,
            colors=[:black], linestyles=[:solid], text_scaling=text_scaling, legend=legend, atm=atm, oce=oce
        )

    elseif plot_unstable_range
        # Plot CFL condition range for atmosphere or ocean
        variable_dict = get_var_dict()
        color_dict, linestyle_dict = get_color_dict()
        physical_values[:delta_t_min] = 100

        delta_z = 10 .^ LinRange(log10(0.1), log10(10), 50)
        delta_t = 10 .^ LinRange(log10(1), log10(100), 50)
        n_zs_atm = atm ? Int.(round.((physical_values[:h_atm] - physical_values[:z_0numA]) ./ reverse(delta_z))) : nothing
        n_ts_atm = atm ? Int.(round.(physical_values[:delta_t_min] ./ reverse(delta_t))) : nothing
        n_zs_oce = oce ? Int.(round.((physical_values[:h_oce] - physical_values[:z_0numO]) ./ reverse(delta_z))) : nothing
        n_ts_oce = oce ? Int.(round.(physical_values[:delta_t_min] ./ reverse(delta_t))) : nothing

        xscale = :log10
        yscale = :log10
        xticks = [1, 10, 100]
        yticks = [0.1, 1, 10]

        a_is = !isempty(a_is) ? a_is : [physical_values[:a_i]]

        plot()
        for a_i in a_is
            physical_values[:a_i] = a_i
            if atm
                unstable_matrix_atm, theoretical_vals_matrix_atm = stability_check(physical_values,
                    n_zs_atm, n_ts_atm, "n_atm", "n_t_atm"
                )
                # Continue here!:) fix this method call. Would be nice to be able to plot the ocean and atmosphere in the same and potentially
                # With the theoretical limit and direction.
                delta_z, _, unstable_matrix_atm, _, _ = handle_variable(
                    n_zs_atm, "n_atm", nothing, unstable_matrix_atm, physical_values;
                    dims=1
                )
                delta_t, _, unstable_matrix_atm, _, _ = handle_variable(
                    n_ts_atm, "n_t_atm", nothing, unstable_matrix_atm, physical_values;
                    dims=2
                )
                theoretical_vals_matrix_atm = reverse(theoretical_vals_matrix_atm, dims=(1, 2))
                plot_delta_z_delta_t(unstable_matrix_atm, theoretical_vals_matrix_atm, delta_z, delta_t, L"$\Delta z^A$", L"$\Delta t^A$",
                    xscale=xscale, yscale=yscale, xticks=xticks, yticks=yticks,
                    color=color_dict[round(a_i, digits=1)], a_i=a_i, legend=legend
                )
            elseif oce
                unstable_matrix_oce, theoretical_vals_matrix_oce = stability_check(physical_values,
                    n_zs_oce, n_ts_oce, "n_oce", "n_t_oce"
                )
                delta_z, theoretical_vals_matrix_oce, unstable_matrix_oce, _, _ = handle_variable(
                    n_zs_oce, "n_oce", unstable_matrix_oce, nothing, physical_values;
                    dims=1
                )
                delta_t, _, unstable_matrix_oce, _, _ = handle_variable(
                    n_ts_oce, "n_t_oce", unstable_matrix_oce, nothing, physical_values;
                    dims=2
                )
                theoretical_vals_matrix_atm = reverse(theoretical_vals_matrix_oce, dims=(1, 2))
                plot_delta_z_delta_t(unstable_matrix_oce, theoretical_vals_matrix_oce, delta_z, delta_t, L"$\Delta z^O$", L"$\Delta t^O$",
                    xscale=xscale, yscale=yscale, xticks=xticks, yticks=yticks,
                    color=color_dict[round(a_i, digits=1)], a_i=a_i, legend=legend
                )
            end
        end
        display(current())
    end

    # conv_facs_analytic = nothing
    # param_analytic = nothing

    # # Plot with respect to a_i
    # var1s = 0:0.1:1

    # # # Plot with respecto to a_i and other parameters
    # # var1s = [0.1, 0.3, 0.5] # a_i
    # # # var2s = [10, 100, 1000, 10000, 100000, 1000000] #t_max
    # # # var2s = [15, 50, 100, 200, 300, 500, 700, 1500] # h_atm
    # # # var2s = [5, 10, 50, 100, 125, 150, 175, 200] # h_oce
    # # # var2s = [20, 200, 2000, 20000] # n_atm
    # # # var2s = [5, 50, 500, 5000, 50000] # n_oce
    # # # var2s = [10, 50, 100, 125, 150, 175, 200] # L_AI
    # # # ai_to_plot = 3
    # # # colors = [:blue, :red, :green]
    # # # var2s = [0.0007, 0.0008, 0.0009, 0.001, 0.0011, 0.0012, 0.0013] # C_AO
    # # # var2s = [0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5] # u_atm and u_oce
    # # # var2s = [260, 262, 264, 266, 268, 270, 273, 276] # T_atm_ini
    # # # var2s = [270, 271, 272, 273, 274, 275, 276] # T_oce_ini
    # # # var2s = [260, 262, 264, 266, 268, 270, 273] # T_ice_ini
    # # var2s = 10 .^ LinRange(log10(1), log10(1000), 10)# n_t_atm
    # # # var2s = 10 .^ LinRange(log10(1), log10(1000), 10)# n_t_oce

    # # xticks = [1, 10, 100, 1000] #var2s
    # # var1_name = "a_i"
    # # var2_name = "n_t_atm"
    # # var2_plot_name = L"$\Delta t^A$"

    # # xscale = :log10
    # # legend = :topright

    # # Run again. See what happens with this better definition of divergence.
    # # Evaluating wrt delta z and delta t. Note that delta z should be first argument.
    # # delta_z = 10 .^ LinRange(log10(0.1), log10(10), 5)
    # # n_atms = Int.(round.((physical_values[:h_atm] - physical_values[:z_0numA]) ./ reverse(delta_z)))
    # # var1s = unique(n_atms)
    # # physical_values[:Δt_min] = 100
    # # delta_t = 10 .^ LinRange(log10(1), log10(100), 5)
    # # n_t_atms = Int.(round.(physical_values[:Δt_min] ./ reverse(delta_t)))
    # # var2s = unique(n_t_atms)

    # # delta_z = 10 .^ LinRange(log10(0.1), log10(10), 10)
    # # n_oces = Int.(round.((physical_values[:h_oce] - physical_values[:z_0numO]) ./ delta_z))
    # # var1s = unique(n_oces)
    # # physical_values[:Δt_min] = 100
    # # delta_t = 10 .^ LinRange(log10(1), log10(100), 10)
    # # n_t_oces = Int.(round.(physical_values[:Δt_min] ./ delta_t))
    # # var2s = unique(n_t_oces)

    # # yticks = [0.1, 1, 10]#var2s
    # # xticks = [1, 10, 100]
    # # var1_name = "n_atm"
    # # var2_name = "n_t_atm"
    # # var1_plot_name = L"$\Delta z^A$"
    # # var2_plot_name = L"$\Delta t^A$"
    # # xscale = :log10
    # # yscale = :log10

    # # # Computing the convergence factor with slightly different commands.
    # # conv_facs_atm, conv_facs_oce, param_analytic, conv_facs_analytic = get_conv_facs_two_variables(physical_values, var1s, var2s, var1_name, var2_name, analytic=true, log_scale=(xscale == :log10))
    # # # conv_facs_atm, conv_facs_oce, param_analytic, conv_facs_analytic, C_AIs, C_AIs_analytic = get_conv_facs_two_variables(physical_values, var1s, var2s, var1_name, var2_name, analytic=true, log_scale=(xscale == :log10))
    # # # conv_facs_atm, conv_facs_oce = get_conv_facs_two_variables(physical_values, var1s, var2s, var1_name, var2_name, analytic=false, log_scale=(xscale == :log10))
    # # conv_facs_atm, conv_facs_oce = get_conv_facs_two_variables(physical_values, var1s, var2s, var1_name, var2_name, analytic=false, a_i=a_i) # for delta z, delta t dependence.
    # # conv_facs_atm, conv_facs_oce = one_iter_update(physical_values, var1s, var2s, var1_name, var2_name)

    # conv_facs_numeric = conv_facs_atm

    # # Special treatment of n_atm and n_oce.
    # # if var1_name == "n_atm"
    # #     var1s = (physical_values[:h_atm] - physical_values[:z_0numA]) ./ reverse(var1s)
    # #     conv_facs_numeric = reverse(conv_facs_numeric, dims=1)
    # # end
    # # if var2_name == "n_atm"
    # #     var2s = (physical_values[:h_atm] - physical_values[:z_0numA]) ./ reverse(var2s)
    # #     conv_facs_numeric = reverse(conv_facs_numeric, dims=2)
    # #     if !isnothing(conv_facs_analytic) && !isnothing(param_analytic)
    # #         param_analytic = reverse((physical_values[:h_atm] - physical_values[:z_0numA]) ./ param_analytic)
    # #         conv_facs_analytic = reverse(conv_facs_analytic, dims=2)
    # #     end
    # #     xticks = reverse((physical_values[:h_atm] - physical_values[:z_0numA]) ./ xticks)
    # # end
    # # if var1_name == "n_oce"
    # #     var1s = (physical_values[:h_oce] - physical_values[:z_0numO]) ./ reverse(var1s)
    # #     conv_facs_numeric = reverse(conv_facs_numeric, dims=1)
    # # end
    # # if var2_name == "n_oce"
    # #     var2s = (physical_values[:h_oce] - physical_values[:z_0numO]) ./ reverse(var2s)
    # #     conv_facs_numeric = reverse(conv_facs_numeric, dims=2)
    # #     if !isnothing(conv_facs_analytic) && !isnothing(param_analytic)
    # #         param_analytic = reverse((physical_values[:h_oce] - physical_values[:z_0numO]) ./ param_analytic)
    # #         conv_facs_analytic = reverse(conv_facs_analytic, dims=2)
    # #     end
    # #     xticks = reverse((physical_values[:h_oce] - physical_values[:z_0numO]) ./ xticks)
    # # end

    # # # Special treatment for n_t_atm and n_t_oce. Only have as second arguments, so only have to handle this.
    # # if var2_name == "n_t_atm" || var2_name == "n_t_oce"
    # #     var2s = physical_values[:Δt_min] ./ reverse(var2s)
    # #     conv_facs_numeric = reverse(conv_facs_numeric, dims=2)
    # #     if var1_name == "a_i"
    # #         xticks = physical_values[:Δt_min] ./ reverse(xticks)
    # #     end
    # #     if !isnothing(conv_facs_analytic) && !isnothing(param_analytic)
    # #         param_analytic = reverse(physical_values[:Δt_min] ./ param_analytic)
    # #         conv_facs_analytic = reverse(conv_facs_analytic, dims=2)
    # #     end
    # # end

    # # # Special treatment for L_AI
    # # if var2_name == "L_AI"
    # #     var2s = sort(C_AIs[ai_to_plot, :])
    # #     param_analytic = sort(C_AIs_analytic[ai_to_plot, :])

    # #     rounded_x = round.(var2s, digits=3)
    # #     rounded_and_unique_x = unique(rounded_x)
    # #     # To store the indices of the first occurrences
    # #     first_indices = []

    # #     # Set to track already seen values
    # #     seen_values = Set()

    # #     # Iterate over the vector
    # #     for (i, value) in enumerate(rounded_x)
    # #         if !(value in seen_values)
    # #             push!(first_indices, i)
    # #             push!(seen_values, value)
    # #         end
    # #     end
    # #     x_ticks = var2s[first_indices]
    # #     xticks = (x_ticks, rounded_and_unique_x)
    # #     sorted_C_AI_indices = sortperm(C_AIs[ai_to_plot, :])
    # #     conv_facs_numeric = conv_facs_numeric[ai_to_plot, :][sorted_C_AI_indices]
    # #     var1s = var1s[ai_to_plot]

    # #     sorted_C_AI_indices_analytic = sortperm(C_AIs_analytic[ai_to_plot, :])
    # #     conv_facs_analytic = conv_facs_analytic[ai_to_plot, :][sorted_C_AI_indices_analytic]
    # # end

    # # Plot wrt ice

    # # # Plot wrt ice and one other param
    # # plot_wrt_a_i_and_one_param(conv_facs_numeric, var1s, var2s, var2_plot_name, conv_facs_analytic=conv_facs_analytic, param_analytic=param_analytic, xscale=xscale, xticks=xticks, scale_text_position=3, text_start_at_bottom=1, legend=legend) # For plotting wrt many params

    # # # Plot wrt one parameter (use for C_AI)
    # # # plot_wrt_one_param(conv_facs_numeric, var2s, var2_plot_name, conv_facs_analytic=conv_facs_analytic, param_analytic=param_analytic, xticks=xticks, log_conv_fac=false, ylim=:auto, color=colors[ai_to_plot]) # For plotting wrt C_AI

    # # # Plot only the numerical convergence factor
    # # # plot_wrt_a_i_and_one_param(conv_facs_numeric, var1s, var2s, var2_plot_name, xscale=xscale, xticks=xticks) # For plotting only the numerical. (for instance with regards to delta t)

    # # Plot wrt delta t and delta z
    # # plot_delta_z_delta_t(conv_facs_numeric, var1s, var2s, var1_plot_name, var2_plot_name, xscale=xscale, yscale=yscale, xticks=xticks, yticks=yticks, color=:blue, a_i=a_i)

    # # # Set appropriate backend for plotting
    # # plotly()
    # # gr()

    # # scatter(x, conv_facs_oce[i, :], xticks=x, xscale=:log10, yformatter=:scientific, label="aᴵ=$ai", markersize=5, xlabel=x_axis_label, ylabel="ρᴼ", legend=:topright, color=colors[i])
    # # plot!(x, conv_facs_oce[i, :], label="", linewidth=2, color=colors[i])

    # # # Plot convergence factor as a function of deltaz delta t quotient
    # # Δzᴼs = reverse((H_O - z_0numO) ./ n_oces) # For delta_z
    # # quotient = Δzᴼs ./ (Δt_min / n_t_oce)
    # # conv_facs_atm[.!isinf.(conv_facs_atm)] .= 0
    # # conv_facs_atm[isinf.(conv_facs_atm)] .= 1
    # # heatmap(conv_facs_atm', xticks=(1:20:length(a_is), round.(a_is[1:20:length(a_is)], digits=1)), yticks=(1:20:length(quotient), round.(quotient[1:20:length(quotient)], digits=4)), color=:viridis,
    # #     xlabel="aᴵ", ylabel="Δzᴼ/Δtᴼ",
    # #     clims=(-1, 1))  # Set color limits    println(conv_facs_atm)
    # # # plot3d(a_is, quotient, conv_facs_atm', label="atm", color=:blue, markersize=5, xlabel="aᵢ", ylabel="", zlabel="ρ", legend=:topright)
    # # # # surface!(a_is[1:length(conv_facs_oce[:, 1])], t_maxs[1:length(conv_facs_oce[1, :])], conv_facs_oce', label="oce", color=:green, markersize=5)
    # # display(current())
end;






