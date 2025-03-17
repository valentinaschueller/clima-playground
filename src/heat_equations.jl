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
            original_file = joinpath(cs.dirs.artifacts, "checkpoint", "checkpoint_" * Interfacer.name(sim) * "_$time.hdf5")
            new_file = joinpath(cs.dirs.artifacts, "checkpoint", "checkpoint_" * Interfacer.name(sim) * "_$iter" * "_$time.hdf5")
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

function solve_coupler!(cs::Interfacer.CoupledSimulation; iterate=10, parallel=false, print_conv=false, plot_conv=false, return_conv=false, analytic_conv_fac_value=nothing, combine=false, atm=true, oce=true, legend=:right)
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
    starting_temp_atm = parent(cs.model_sims.atmos_sim.Y_init)
    starting_temp_oce = parent(cs.model_sims.ocean_sim.Y_init)
    starting_temp_ice = parent(cs.model_sims.ice_sim.Y_init)
    upper_limit_temp = maximum([maximum(starting_temp_oce), maximum(starting_temp_atm), maximum(starting_temp_ice)])
    lower_limit_temp = minimum([minimum(starting_temp_oce), minimum(starting_temp_atm), minimum(starting_temp_ice)])

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
function coupled_heat_equations(; iterate=10, parallel=false, boundary_mapping="mean", values=Dict{Symbol,Int}(), print_conv_facs_iter=false, plot_conv_facs_iter=false, analytic_conv_fac=false, atm=true, oce=true, combine=false, plot_unstable_range=false, a_is=[], var_name=nothing, xscale=:identity, legend=:right, log_conv_fac=false, xticks=nothing, text_scaling=(1, 5))
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

        delta_z = 10 .^ LinRange(log10(0.001), log10(1), 50)
        delta_t = 10 .^ LinRange(log10(1), log10(100), 50)
        n_zs_atm = atm ? Int.(round.((physical_values[:h_atm] - physical_values[:z_0numA]) ./ reverse(delta_z))) : nothing
        n_ts_atm = atm ? Int.(round.(physical_values[:delta_t_min] ./ reverse(delta_t))) : nothing
        n_zs_oce = oce ? Int.(round.((physical_values[:h_oce] - physical_values[:z_0numO]) ./ reverse(delta_z))) : nothing
        n_ts_oce = oce ? Int.(round.(physical_values[:delta_t_min] ./ reverse(delta_t))) : nothing

        xscale = :log10
        yscale = :log10
        xticks = [1, 10, 100]
        yticks = [0.001, 0.01, 0.1]

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
                delta_z, unstable_matrix_oce, _, _, _ = handle_variable(
                    n_zs_oce, "n_oce", unstable_matrix_oce, nothing, physical_values;
                    dims=1
                )
                delta_t, unstable_matrix_oce, _, _, _ = handle_variable(
                    n_ts_oce, "n_t_oce", unstable_matrix_oce, nothing, physical_values;
                    dims=2
                )
                theoretical_vals_matrix_oce = reverse(theoretical_vals_matrix_oce, dims=(1, 2))
                plot_delta_z_delta_t(unstable_matrix_oce, theoretical_vals_matrix_oce, delta_z, delta_t, L"$\Delta z^O$", L"$\Delta t^O$",
                    xscale=xscale, yscale=yscale, xticks=xticks, yticks=yticks,
                    color=color_dict[round(a_i, digits=1)], a_i=a_i, legend=legend
                )
            end
        end
        display(current())
    end
end;






