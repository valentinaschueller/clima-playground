"""
Returns the ylabel.

**Arguments:**

-`plot_numeric::Boolean`: Whether to plot the numerical convergence factor.
-`conv_fac_analytic::Array`: Containing the analytical convergence factors for a variable.

"""
function get_ylabel_or_nothing_to_plot(plot_numeric, conv_facs_analytic)
    if plot_numeric && !isnothing(conv_facs_analytic)
        ylabel = L"$\hat{\rho}$, $\rho$"
    elseif plot_numeric && isnothing(conv_facs_analytic)
        ylabel = L"$\rho$"
    elseif !plot_numeric && !isnothing(conv_facs_analytic)
        ylabel = L"$\hat{\rho}$"
    elseif !plot_numeric && isnothing(conv_facs_analytic)
        return nothing
    end
    return ylabel
end

"""
Plots the convergence factors.

**Arguments:**

-`conv_facs_oce::Array`: The ocean convergence factors.
-`conv_facs_atm::Array`: The atmosphere convergence factors.
-`param_numeric::Array`: The values for the variable.
-`param_name::String`: The name of the variable.

**Optional Keyword Arguments:**

-`conv_facs_analytic::Array`: The analytic convergence factors, default: `nothing`.
-`param_analytic::Array`: The analytical values for the variable, default: `nothing`.
-`xscale::Symbol`: Plotting argument for x-scale, default: `:identity`.
-`yscale::Symbol`: Plotting argument for y-scale, default: `:identity`.
-`label::String`: Label for the legend, for instance `a_i=0.4`, default: `""`.
-`color::Symbol`: The color of the plot lines or markers, default: `:black`.
-`linestyle::Symbol`: Style of the plot lines, default: `:solid`.
-`xticks::Symbol or Array`: Defines the tick marks on the x-axis, default: `:auto`.
-`yticks::Symbol or Array`: Defines the tick marks on the y-axis, default: `:auto`.
-`ylim::Symbol or Tuple`: The limits for the y-axis, default: `:auto`.
-`legend::Symbol`: Legend position, default: `:right`.
-`compute_atm_conv_fac::Boolean`: Whether to plot the atmosphere convergence factor, default: `true`.
-`compute_oce_conv_fac::Boolean`: Whether to plot the ocean convergence factor, default: `true`.

"""
function plot_data(
    conv_facs_oce,
    conv_facs_atm,
    param_numeric,
    param_name;
    conv_facs_analytic=nothing,
    param_analytic=nothing,
    xscale=:identity,
    yscale=:identity,
    label="",
    color=:black,
    linestyle=:solid,
    xticks=:auto,
    yticks=:auto,
    ylim=:auto,
    legend=:right,
    compute_atm_conv_fac=true,
    compute_oce_conv_fac=true,
)
    ylabel = get_ylabel_or_nothing_to_plot(
        !isnothing(conv_facs_oce) || !isnothing(conv_facs_atm),
        conv_facs_analytic,
    )
    if isnothing(ylabel)
        return nothing
    end
    if !isnothing(conv_facs_analytic)
        param_analytic_new =
            yscale == :log10 ? param_analytic[conv_facs_analytic.>0] : param_analytic
        conv_facs_analytic_new =
            yscale == :log10 ? conv_facs_analytic[conv_facs_analytic.>0] :
            conv_facs_analytic
        plot!(
            param_analytic_new,
            conv_facs_analytic_new,
            label="analytical" * label,
            color=color,
            xscale=xscale,
            yscale=yscale,
            linestyle=linestyle,
            linewidth=2,
            legend=legend,
            xticks=xticks,
            yticks=yticks,
            ylim=ylim,
            legend_background_color=RGBA(1, 1, 1, 0.5),
            legendfontsize=12,
        )
    end
    if compute_oce_conv_fac
        plot!(
            param_numeric,
            conv_facs_oce,
            label="numerical (oce)" * label,
            color=color,
            xscale=xscale,
            yscale=yscale,
            linestyle=linestyle,
            markershape=:circle,
            linewidth=2,
            legend=legend,
            xticks=xticks,
            yticks=yticks,
            ylim=ylim,
            legend_background_color=RGBA(1, 1, 1, 0.5),
            legendfontsize=12,
        )
    end
    if compute_atm_conv_fac
        plot!(
            param_numeric,
            conv_facs_atm,
            label="numerical (atm)" * label,
            color=color,
            xscale=xscale,
            yscale=yscale,
            linestyle=linestyle,
            markershape=:x,
            linewidth=2,
            legend=legend,
            xticks=xticks,
            yticks=yticks,
            ylim=ylim,
            legend_background_color=RGBA(1, 1, 1, 0.5),
            legendfontsize=12,
        )
    end
    xlabel!(param_name)
    ylabel!(ylabel)
end

"""
Computes and plots the convergence factor for different a_i and wrt one parameter. Some special treatment
for instabilities has been added.

**Arguments:**

-`conv_facs_oce::Array`: The ocean convergence factors.
-`conv_facs_atm::Array`: The atmosphere convergence factors.
-`param_numeric::Array`: The values for the variable.
-`param_name::String`: The name of the variable.

**Optional Keyword Arguments:**

-`conv_facs_analytic::Array`: The analytic convergence factors, default: `nothing`.
-`param_analytic::Array`: The analytical values for the variable, default: `nothing`.
-`xscale::Symbol`: Plotting argument for x-scale, default: `:identity`.
-`yscale::Symbol`: Plotting argument for y-scale, default: `:identity`.
-`colors::Array`: The colors of the plot lines and markers, default: `[:blue, :red, :green]`.
-`linestyles::Array`: Styles of the plot lines, default: `[:solid, :dash, :dot]`
-`xticks::Symbol or Array`: Defines the tick marks on the x-axis, default: `:auto`.
-`yticks::Symbol or Array`: Defines the tick marks on the y-axis, default: `:auto`.
-`text_scaling::Tuple`: If the model goes unstable, there is plot in the text, this allows for different text positions, default: `(1, 5)`.
-`legend::Symbol`: Legend position, default: `:right`.
-`compute_atm_conv_fac::Boolean`: Whether to plot the atmosphere convergence factor, default: `true`.
-`compute_oce_conv_fac::Boolean`: Whether to plot the ocean convergence factor, default: `true`.

"""
function plot_wrt_a_i_and_one_param(
    conv_facs_oce,
    conv_facs_atm,
    a_is,
    param_numeric,
    param_name;
    conv_facs_analytic=nothing,
    param_analytic=nothing,
    xscale=:identity,
    yscale=:identity,
    colors=[:blue, :red, :green],
    linestyles=[:solid, :dash, :dot],
    xticks=:auto,
    yticks=:auto,
    text_scaling=(1, 5),
    legend=:right,
    compute_atm_conv_fac=true,
    compute_oce_conv_fac=true,
)
    gr()
    plot()
    # If there is instabilities, handle the text
    unstable_oce_indices = isinf.(conv_facs_oce)
    conv_facs_oce[unstable_oce_indices] .= NaN

    unstable_atm_indices = isinf.(conv_facs_atm)
    conv_facs_atm[unstable_atm_indices] .= NaN

    analytic = !isnothing(conv_facs_analytic)

    lowerbound, upperbound =
        get_bounds(yscale, conv_facs_analytic, conv_facs_oce, conv_facs_atm)

    for (i, ai) in enumerate(a_is)
        conv_facs_oce_i = conv_facs_oce[i, :]
        conv_facs_atm_i = conv_facs_atm[i, :]
        conv_facs_analytic_i = analytic ? conv_facs_analytic[i, :] : nothing

        x_text_oce = param_numeric[unstable_oce_indices[i, :]]
        x_text_atm = param_numeric[unstable_atm_indices[i, :]]
        param_numeric_i = copy(param_numeric)
        param_numeric_i[unstable_oce_indices[i, :]] .= NaN
        param_numeric_i[unstable_atm_indices[i, :]] .= NaN

        if yscale == :log10
            y_text_oce = [
                10^(
                    log10(lowerbound) + (
                        (i - text_scaling[1]) * (log10(upperbound) - log10(lowerbound)) / text_scaling[2]
                    )
                ) for _ in x_text_oce
            ]
            y_text_atm = [
                10^(
                    log10(lowerbound) + (
                        (i - text_scaling[1]) * (log10(upperbound) - log10(lowerbound)) / text_scaling[2]
                    )
                ) for _ in setdiff(x_text_atm, x_text_oce)
            ]
        else
            y_text_oce = [
                lowerbound +
                ((i - text_scaling[1]) * (upperbound - lowerbound) / text_scaling[2])
                for _ in x_text_oce
            ]
            y_text_atm = [
                lowerbound +
                ((i - text_scaling[1]) * (upperbound - lowerbound) / text_scaling[2])
                for _ in setdiff(x_text_atm, x_text_oce)
            ]
        end

        plot_data(
            conv_facs_oce_i,
            conv_facs_atm_i,
            param_numeric_i,
            param_name,
            conv_facs_analytic=conv_facs_analytic_i,
            param_analytic=param_analytic,
            xscale=xscale,
            yscale=yscale,
            label=(param_name !== L"$a^I$") ? ", aᴵ = " * "$ai" : "",
            color=colors[i],
            linestyle=linestyles[i],
            xticks=xticks,
            yticks=yticks,
            legend=legend,
            compute_atm_conv_fac=compute_atm_conv_fac,
            compute_oce_conv_fac=compute_oce_conv_fac,
        )
        for (k, txt) in enumerate(y_text_oce)
            annotate!(
                x_text_oce[k],
                txt,
                text("oce \u26A1", 12, colors[i], :left, rotation=90),
            )
        end
        for (k, txt) in enumerate(y_text_atm)
            annotate!(
                x_text_atm[k],
                txt,
                text("atm \u26A1", 12, colors[i], :left, rotation=90),
            )
        end
    end
    display(current())
end

"""
Computes the upper and lower bounds for the unstable text.

**Arguments:**

-`yscale::Symbol`: Plotting argument for y-scale.
-`conv_facs_analytic::Array`: The analytic convergence factors.
-`finite_oce_vals::Array`: The finite values of the ocean convergence factors.
-`finite_atm_vals::Array`: The finite values of the atmosphere convergence factors.

"""
function get_bounds(yscale, conv_facs_analytic, finite_oce_vals, finite_atm_vals)
    matrices = filter(!isnothing, [conv_facs_analytic, finite_oce_vals, finite_atm_vals])
    filtered_matrices = [filter(!isnan, matrix) for matrix in matrices]
    filtered_matrices =
        yscale == :log10 && !isempty(filtered_matrices) ?
        [filter(x -> x != 0, matrix) for matrix in filtered_matrices] : filtered_matrices
    upperbound =
        !isempty(filtered_matrices) ?
        maximum(maximum.(filter(!isempty, filtered_matrices))) : 1
    lowerbound =
        !isempty(filtered_matrices) ?
        minimum(minimum.(filter(!isempty, filtered_matrices))) : 0
    return lowerbound, upperbound
end

"""
Plots the unstable model regime.

**Arguments:**

-`unstable_matrix: Array`: Containing the indices (with respect to Δz and Δt) at which the model becomes unstable.
-`Δzs::Array`: Range to vary Δzᴬ or Δzᴼ over.
-`Δts::Array`: Range to vary Δtᴬ or Δtᴼ over.
-`Δz_name::String`: Either Δzᴬ or Δzᴼ.
-`Δt_name::String`: Either Δtᴬ or Δtᴼ.

**Optional Keyword Arguments:**

-`xscale::Symbol`: Plotting argument for x-scale, default: `:identity`.
-`yscale::Symbol`: Plotting argument for y-scale, default: `:identity`.
-`xticks::Symbol or Array`: Defines the tick marks on the x-axis, default: `:auto`.
-`yticks::Symbol or Array`: Defines the tick marks on the y-axis, default: `:auto`.
-`a_i::Float64`: Sea ice concentration, default: `Float64(0.5)`.
-`color::Symbol`: Color for the plotted regime, default: `:green`.
-`legend::Symbol`: Legend position, default: `:right`.

"""
function plot_Δz_Δt(
    unstable_matrix,
    Δzs,
    Δts,
    Δz_name,
    Δt_name;
    xscale=:identity,
    yscale=:identity,
    xticks=:auto,
    yticks=:auto,
    a_i=Float64(0.5),
    color=:green,
    legend=:right,
)
    t_values = []
    pre_index = 0
    for (i, Δz) in enumerate(Δzs)
        first_inf_index = findfirst(isinf, unstable_matrix[i, :])

        if isnothing(first_inf_index)
            push!(t_values, NaN)
        else
            first_t_inf = Δts[first_inf_index]
            push!(t_values, first_t_inf)
            if first_inf_index == pre_index
                t_values[i-1] = NaN
            end
        end
        pre_index = first_inf_index
    end

    not_nan_indices = .!isnan.(t_values)
    if not_nan_indices[end]
        last_t = t_values[end]
        index_for_last_t = findall(x -> x == last_t, Δts)[1]
        remaining_t = Δts[index_for_last_t:end]
        t_values = vcat(t_values[not_nan_indices], remaining_t)
        Δzs_new = vcat(Δzs[not_nan_indices], (yticks[end] + 1) * ones(length(remaining_t)))
    else
        t_values = t_values[not_nan_indices]
        Δzs_new = Δzs[not_nan_indices]
    end
    plot!(
        t_values,
        Δzs_new,
        linewidth=2,
        xscale=xscale,
        yscale=yscale,
        xticks=xticks,
        yticks=yticks,
        xlabel=Δt_name,
        ylabel=Δz_name,
        label="Unstable regime, aᴵ=$a_i",
        color=color,
        fillrange=minimum(Δzs),
        fillalpha=0.3,
        xlim=(xticks[1], xticks[end]),
        ylim=(yticks[1], yticks[end]),
        legend=legend,
    )

    scatter!(t_values, Δzs_new, markershape=:circle, label="", color=color)
    if Δz_name == L"$\Delta z^A$"
        t1 = 2:1:7
        t2 = 20:10:70
        plot!(
            t1,
            (t1 .^ (1 / 2)) .* 10^-2.5,
            color=:black,
            label="",
            xscale=xscale,
            yscale=yscale,
        )
        plot!(t2, t2 ./ 1000, color=:black, label="", xscale=xscale, yscale=yscale)
    else
        t2 = 10 .^ LinRange(log10(5), log10(20), 50)
        plot!(
            t2,
            t2 .* 10 .^ -2.7,
            color=:black,
            label="",
            xscale=xscale,
            yscale=yscale,
        )
    end


    if Δz_name == L"$\Delta z^A$"
        annotate!(
            4,
            (5 .^ (1 / 2)) .* 10^(-2.5),
            Plots.text(L"slope = 1/2", 10, :black, rotation=18),
        )
        annotate!(40, 50 ./ 1000, Plots.text(L"slope = 1", 10, :black, rotation=34))
    else
        annotate!(
            t2[25],
            t2[25] .* 10^-2.6,
            Plots.text(L"slope = 1", 10, :black, rotation=34),
        )
    end

end

function plot_ρ_over_k(iterations; parallel=false, analytic_conv_fac=true, combine_ρ_parallel=false, compute_atm_conv_fac=true, compute_oce_conv_fac=true, legend=:right)
    cs, conv_fac_atm, conv_fac_oce = coupled_heat_equations(iterations=iterations, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000))
    physical_values = cs.model_sims.atmos_sim.params

    if analytic_conv_fac
        analytic_conv_fac_value = compute_ρ_analytical(physical_values)
        combine_ρ_parallel = true
    end

    if parallel && combine_ρ_parallel
        conv_fac_atm, conv_fac_oce = update_ρ_parallel(conv_fac_atm, conv_fac_oce)
        ylabel = L"$\rho_{k+1}\times\rho_k$"
    else
        ylabel = L"$\rho_k$"
    end

    if analytic_conv_fac
        ylabel *= L", $\hat{\rho}$"
    end

    gr()
    plot()
    color_dict, _ = get_color_dict()
    color = color_dict[round(cs.model_sims.atmos_sim.params.a_i, digits=1)]
    k_atm = 2:length(conv_fac_atm)+1
    k_oce = 2:length(conv_fac_oce)+1
    if compute_atm_conv_fac
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
    end
    if compute_oce_conv_fac
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
            ylim=(
                0,
                maximum([
                    maximum(conv_fac_atm[.!isnan.(conv_fac_atm)]),
                    maximum(conv_fac_oce[.!isnan.(conv_fac_oce)]),
                ]) * 1.2,
            ),
        )
    end
    if analytic_conv_fac
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

function plot_unstable_range(component; a_is=[])
    physical_values = define_realistic_vals()
    correct_for_a_i!(physical_values)
    compute_C_AO!(physical_values)
    physical_values[:Δt_min] = 100
    color_dict, _ = get_color_dict()

    Δz = 10 .^ LinRange(log10(0.001), log10(1), 50)
    Δt = 10 .^ LinRange(log10(1), log10(100), 50)
    if component == "atm"
        n_zs_atm = Int.(round.((physical_values[:h_atm] - physical_values[:z_0numA]) ./ reverse(Δz)))
        n_ts_atm = Int.(round.(physical_values[:Δt_min] ./ reverse(Δt)))
    elseif component == "oce"
        n_zs_oce = Int.(round.((physical_values[:h_oce] - physical_values[:z_0numO]) ./ reverse(Δz)))
        n_ts_oce = Int.(round.(physical_values[:Δt_min] ./ reverse(Δt)))
    else
        error("Component must be 'atm' or 'oce'.")
    end

    xscale = :log10
    yscale = :log10
    xticks = [1, 10, 100]
    yticks = [0.001, 0.01, 0.1]
    legend = :right

    a_is = !isempty(a_is) ? a_is : [physical_values[:a_i]]

    plot()
    for a_i in a_is
        physical_values[:a_i] = a_i
        if component == "atm"
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
        else
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


function plot_ρ_over_a_i(iterations=10, analytic_conv_fac=true)
    xscale = :identity
    yscale = :log10
    legend = :right
    yticks = :auto
    a_is = 0:0.1:1
    text_scaling = (1, 5)
    physical_values = define_realistic_vals()
    params = Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000)
    merge!(physical_values, params)
    conv_facs_atm, conv_facs_oce, param_analytic, conv_facs_analytic =
        get_conv_facs_one_variable(
            physical_values,
            a_is,
            "a_i",
            iterations=iterations,
            analytic=analytic_conv_fac,
            log_scale=(xscale == :log10),
        )
    xticks = a_is
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
        compute_atm_conv_fac=true,
        compute_oce_conv_fac=false,
    )
end

function plot_ρ_over_var(iterations, var_name; a_is=[], analytic_conv_fac=true, xticks=nothing, xscale=:identity)
    yscale = :log10
    legend = :right
    yticks = :auto
    text_scaling = (1, 5)
    physical_values = define_realistic_vals()
    params = Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000)
    merge!(physical_values, params)
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
            compute_atm_conv_fac=true,
            compute_oce_conv_fac=false,
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
            compute_atm_conv_fac=true,
            compute_oce_conv_fac=false,
        )
    end
end
