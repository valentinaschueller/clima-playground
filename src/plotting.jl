function get_ylabel_or_nothing_to_plot(log_conv_fac, plot_numeric, conv_facs_analytic)
    if log_conv_fac
        if plot_numeric && !isnothing(conv_facs_analytic)
            ylabel = L"$log(̂\hat{\rho})$, $log(\rho)$"
        elseif plot_numeric && isnothing(conv_facs_analytic)
            ylabel = L"$log(\rho)$"
        elseif !plot_numeric && !isnothing(conv_facs_analytic)
            ylabel = L"$log(̂\hat{\rho})$"
        elseif !plot_numeric && isnothing(conv_facs_analytic)
            return nothing
        end
    else
        if plot_numeric && !isnothing(conv_facs_analytic)
            ylabel = L"$\hat{\rho}$, $\rho$"
        elseif plot_numeric && isnothing(conv_facs_analytic)
            ylabel = L"$\rho$"
        elseif !plot_numeric && !isnothing(conv_facs_analytic)
            ylabel = L"$\hat{\rho}$"
        elseif !plot_numeric && isnothing(conv_facs_analytic)
            return nothing
        end
    end
    return ylabel
end
function plot_data(conv_facs_oce, conv_facs_atm, param_numeric, param_name; conv_facs_analytic=nothing, param_analytic=nothing, xscale=:identity, log_conv_fac=false, label="", color=:black, linestyle=:solid, xticks=:auto, ylim=:auto, legend=:right, atm=true, oce=true)
    ylabel = get_ylabel_or_nothing_to_plot(log_conv_fac, !isnothing(conv_facs_oce) || !isnothing(conv_facs_atm), conv_facs_analytic)
    if ylabel == nothing
        return nothing
    end
    if !isnothing(conv_facs_analytic)
        # Plot analytic data and a legend
        plot!(param_analytic, conv_facs_analytic, label="analytical" * label, color=color, xscale=xscale, linestyle=linestyle, linewidth=2, legend=legend, xticks=xticks, ylim=ylim, legend_background_color=RGBA(1,1,1,0.5))
    end
    if oce
        plot!(param_numeric, conv_facs_oce, label="numerical (oce)" * label, color=color, xscale=xscale, linestyle=linestyle, markershape=:circle, linewidth=2, legend=legend, ylim=ylim, legend_background_color=RGBA(1,1,1,0.5))
    end
    if atm
        plot!(param_numeric, conv_facs_atm, label="numerical (atm)" * label, color=color, xscale=xscale, linestyle=linestyle, markershape=:x, linewidth=2, legend=legend, ylim=ylim, legend_background_color=RGBA(1,1,1,0.5))
    end
    xlabel!(param_name)
    ylabel!(ylabel)
end

function plot_wrt_a_i_and_one_param(conv_facs_oce, conv_facs_atm, a_is, param_numeric, param_name; conv_facs_analytic=nothing, param_analytic=nothing, xscale=:identity, xticks=xticks, colors=[:blue, :red, :green], linestyles=[:solid, :dash, :dot], text_scaling=(1, 5), legend=:right, log_conv_fac=false, atm=true, oce=true)
    """
    Plot conv factor for different a_i and wrt a parameter. Some special treatment
    for divergence has been added
    """
    gr()
    plot()
    # If there is divergence, to handle the text
    ok_oce_indices = .!isinf.(conv_facs_oce)
    unstable_oce_indices = isinf.(conv_facs_oce)
    conv_facs_oce[unstable_oce_indices] .= NaN

    ok_atm_indices = .!isinf.(conv_facs_atm)
    unstable_atm_indices = isinf.(conv_facs_atm)
    conv_facs_atm[unstable_atm_indices] .= NaN

    analytic = !isnothing(conv_facs_analytic)

    if log_conv_fac
        conv_facs_oce = log.(conv_facs_oce)
        conv_facs_atm = log.(conv_facs_atm)
        conv_facs_analytic = analytic ? log.(conv_facs_analytic) : nothing
    end

    lowerbound, upperbound = get_bounds(conv_facs_analytic, conv_facs_oce, conv_facs_atm, log_conv_fac)

    for (i, ai) in enumerate(a_is)
        conv_facs_oce_i = conv_facs_oce[i, :]
        conv_facs_atm_i = conv_facs_atm[i, :]
        conv_facs_analytic_i = analytic ? conv_facs_analytic[i, :] : nothing

        x_text_oce = param_numeric[unstable_oce_indices[i, :]]
        x_text_atm = param_numeric[unstable_atm_indices[i, :]]
        param_numeric_i = copy(param_numeric)
        param_numeric_i[unstable_oce_indices[i, :]] .= NaN
        param_numeric_i[unstable_atm_indices[i, :]] .= NaN

        y_text_oce = [lowerbound + ((i - text_scaling[1]) * (upperbound - lowerbound) / text_scaling[2]) for _ in x_text_oce]
        y_text_atm = [lowerbound + ((i - text_scaling[1]) * (upperbound - lowerbound) / text_scaling[2]) for _ in setdiff(x_text_atm, x_text_oce)]

        plot_data(conv_facs_oce_i, conv_facs_atm_i, param_numeric_i, param_name, conv_facs_analytic=conv_facs_analytic_i, param_analytic=param_analytic, xscale=xscale, log_conv_fac=log_conv_fac, label=(param_name !== L"$a^I$") ? ", aᴵ = " * "$ai" : "", color=colors[i], linestyle=linestyles[i], xticks=xticks, legend=legend, atm=atm, oce=oce)
        for (k, txt) in enumerate(y_text_oce)
            annotate!(x_text_oce[k], txt, text("oce \u26A1", 12, colors[i], :left, rotation=90))
        end
        for (k, txt) in enumerate(y_text_atm)
            annotate!(x_text_atm[k], txt, text("atm \u26A1", 12, colors[i], :left, rotation=90))
        end
    end
    display(current())
end

function get_bounds(conv_facs_analytic, finite_oce_vals, finite_atm_vals, log_conv_fac)
    matrices=filter(!isnothing, [conv_facs_analytic, finite_oce_vals, finite_atm_vals])
    filtered_matrices = [filter(!isnan, matrix) for matrix in matrices]
    upperbound= !isempty(filtered_matrices) ? maximum(maximum.(filtered_matrices)) : 20
    lowerbound= !isempty(filtered_matrices) ? minimum(minimum.(filtered_matrices)) : -20
    return lowerbound, upperbound
end

function plot_delta_z_delta_t(unstable_matrix, theoretical_vals_matrix, delta_zs, delta_ts, delta_z_name, delta_t_name; xscale=:identity, yscale=:identity, xticks=:auto, yticks=:auto, a_i=0.5, color=:green, legend=:right)
    # Plot convergence factor as a function of deltaz and delta t
    t_values = []
    t_values_theoretical = []
    pre_index = 0
    pre_index_theoretical = 0
    for (i, delta_z) in enumerate(delta_zs)
        first_inf_index = findfirst(isinf, unstable_matrix[i, :])
        first_unstable_index = findfirst(x -> maximum(x) > 1, theoretical_vals_matrix[i, :])

        if isnothing(first_inf_index)
            push!(t_values, NaN)
        else
            first_t_inf = delta_ts[first_inf_index]
            push!(t_values, first_t_inf)
            if first_inf_index == pre_index
                t_values[i-1] = NaN
            end
        end
        if isnothing(first_unstable_index)
            push!(t_values_theoretical, NaN)
        else
            first_t_theoretical = delta_ts[first_unstable_index]
            push!(t_values_theoretical, first_t_theoretical)
            if first_unstable_index == pre_index_theoretical
                t_values_theoretical[i-1] = NaN
            end
        end
        pre_index = first_inf_index
        pre_index_theoretical = first_unstable_index
    end

    not_nan_indices = .!isnan.(t_values)
    if not_nan_indices[end]
        last_t = t_values[end]
        index_for_last_t = findall(x -> x == last_t, delta_ts)[1]
        remaining_t = delta_ts[index_for_last_t:end]
        t_values = vcat(t_values[not_nan_indices], remaining_t)
        delta_zs_new = vcat(delta_zs[not_nan_indices], (yticks[end] + 1) * ones(length(remaining_t)))
    else
        t_values = t_values[not_nan_indices]
        delta_zs_new = delta_zs[not_nan_indices]
    end
    plot!(t_values, delta_zs_new, linewidth=2,
        xscale=xscale, yscale=yscale, xticks=xticks, yticks=yticks, xlabel=delta_t_name,
        ylabel=delta_z_name, label="Unstable regime, aᴵ=$a_i", color=color, fillrange=minimum(delta_zs),
        fillalpha=0.3, xlim=(xticks[1], xticks[end]), ylim=(yticks[1], yticks[end]), legend=legend
    )

    scatter!(t_values, delta_zs_new, markershape=:circle, label="", color=color)
    if delta_z_name == L"$\Delta z^A$"
        t1 = 2:1:7
        t2 = 20:10:70
        plot!(t1, (t1 .^ (1 / 2)) .* 10^-2.5, color=:black, label="", xscale=xscale, yscale=yscale)
        plot!(t2, t2 ./ 1000, color=:black, label="", xscale=xscale, yscale=yscale)
    else
        t2 = 10 .^ LinRange(log10(5), log10(20), 50)
        plot!(t2, t2 .* 10 .^ -2.7, color=:black, label="", xscale=xscale, yscale=yscale)
    end


    not_nan_indices_theoretical = .!isnan.(t_values_theoretical)
    if delta_z_name == L"$\Delta z^A$"
        annotate!(4, (5 .^ (1 / 2)) .* 10^(-2.5), Plots.text(L"slope = 1/2", 10, :black, rotation=18))
        annotate!(40, 50 ./ 1000, Plots.text(L"slope = 1", 10, :black, rotation=34))
    else
        annotate!(t2[25], t2[25] .* 10^-2.6, Plots.text(L"slope = 1", 10, :black, rotation=34))
    end

end