using clima_playground
using Plots
using LaTeXStrings
import ClimaCoupler: Interfacer

function plot_C_H_AO_dependence()
    params = SimulationParameters()
    L_AOs = vec(-200:200)
    C_H_AO = zeros(length(L_AOs))
    for (j, L_AO) in enumerate(L_AOs)
        C_H_AO[j] = compute_C_H_AO(params; L_AO=L_AO)
    end
    plot(
        L_AOs,
        C_H_AO,
        xlabel=L"L_{AO}",
        ylabel=L"C_{H,AO}",
        label="",
        color=:black,
    )
end

function plot_C_H_AI_dependence()
    a_Is = 0:0.01:1
    L_AIs = vec(-200:200)
    params = SimulationParameters()
    C_H_AI = zeros(length(a_Is), length(L_AIs))
    for (i, a_I) in enumerate(a_Is)
        params.a_I = a_I
        for (j, L_AI) in enumerate(L_AIs)
            C_H_AI[i, j] = compute_C_H_AI(params; L_AI=L_AI)
        end
    end
    plot()
    surface(L_AIs, a_Is, C_H_AI, xlabel=L"L_{AI}", ylabel=L"a_I", zlabel=L"C_{H,AI}")
end


function analytical_convergence_factor_dependence()
    νs = range(0, stop=10, length=50)
    ωs = range(0.001, stop=10, length=50)
    a_Is = range(0.001, stop=1, length=10)
    ϱs = zeros(length(νs), length(ωs), length(a_Is))
    params = SimulationParameters()

    # Loop over values
    for (k, a_I) in enumerate(a_Is)
        for (i, ν) in enumerate(νs)
            for (j, ω) in enumerate(ωs)
                params.a_I = a_I
                ϱs[i, j, k] = compute_ϱ_analytical(params; s=ν + im * ω)
            end
        end
        i_max, j_max = Tuple(CartesianIndices(ϱs)[argmax(ϱs)])
        ν_max = νs[i_max]
        ω_max = ωs[j_max]
        println("supremum at ν=$ν_max and ω=$ω_max")
    end

    surface(
        ωs,
        νs,
        ϱs[:, :, 1],
        color=:viridis,
        xlabel="ω",
        ylabel="ν",
        zlabel="ϱ(ν+iω)",
    )
end


function plot_data(
    ϱs_oce,
    ϱs_atm,
    param_numeric,
    param_name;
    ϱs_analytic=nothing,
    param_analytic=nothing,
    kwargs...
)
    ylabel = "ϱ"
    append_label = get(kwargs, :append_label, "")
    if !isnothing(ϱs_analytic)
        plot!(
            param_analytic[ϱs_analytic.>0],
            ϱs_analytic[ϱs_analytic.>0],
            label=L"$ϱ_\mathrm{ana}$" * append_label,
            linewidth=2,
            legend_background_color=RGBA(1, 1, 1, 0.5),
            legendfontsize=12;
            kwargs...
        )
    end
    plot!(
        param_numeric,
        ϱs_atm,
        label=L"$ϱ_\mathrm{num}$" * append_label,
        markershape=:x,
        linewidth=2,
        legend_background_color=RGBA(1, 1, 1, 0.5),
        legendfontsize=12;
        kwargs...
    )
    xlabel!(param_name)
    ylabel!(ylabel)
end


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
    a_I=Float64(0.5),
    color=:green,
    legend=:right,
    kwargs...
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
        Δzs_new;
        linewidth=2,
        xscale=xscale,
        yscale=yscale,
        xticks=xticks,
        yticks=yticks,
        xlabel=Δt_name,
        ylabel=Δz_name,
        label="Unstable regime, aᴵ=$a_I",
        color=color,
        fillrange=minimum(Δzs),
        xlim=(xticks[1], xticks[end]),
        ylim=(yticks[1], yticks[end]),
        legend=legend,
        kwargs...
    )

    scatter!(t_values, Δzs_new, markershape=:circle, label="", color=color)
    if Δz_name == L"$\Delta z^A$"
        t1 = 2:1:7
        t2 = 20:10:70
        plot!(
            t1,
            (t1 .^ (1 / 2)) .* 10^-2.5;
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
            t2 .* 10 .^ -2.7;
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

function plot_unstable_range(component; a_Is=[0.0])
    p = SimulationParameters(Δt_min=100, T_Ib=270.0)

    Δz = 10 .^ LinRange(log10(0.001), log10(1), 50)
    Δt = 10 .^ LinRange(log10(1), log10(100), 50)
    if component == "atm"
        n_zs = Int.(round.((p.h_A - p.z_A0) ./ reverse(Δz)))
        n_ts = Int.(round.(p.Δt_min ./ reverse(Δt)))
        space_field = :n_A
        time_field = :n_t_A
        xlabel = L"$\Delta z^A$"
        ylabel = L"$\Delta t^A$"
    elseif component == "oce"
        n_zs = Int.(round.((p.h_O - p.z_O0) ./ reverse(Δz)))
        n_ts = Int.(round.(p.Δt_min ./ reverse(Δt)))
        space_field = :n_O
        time_field = :n_t_O
        xlabel = L"$\Delta z^O$"
        ylabel = L"$\Delta t^O$"
    else
        error("Component must be 'atm' or 'oce'.")
    end

    xticks = [1, 10, 100]
    yticks = [0.001, 0.01, 0.1]

    plot()
    for a_I in a_Is
        p.a_I = a_I
        p.C_H_AO = compute_C_H_AO(p)
        p.C_H_AI = compute_C_H_AI(p)
        restore_physical_values!(p)
        unstable_matrix = stability_check(p, n_zs, n_ts, space_field, time_field)
        Δz, _, unstable_matrix, _, _ = treat_grid_sizes(
            n_zs,
            space_field,
            nothing,
            unstable_matrix,
            p;
            dims=1,
        )
        Δt, _, unstable_matrix, _, _ = treat_grid_sizes(
            n_ts,
            time_field,
            nothing,
            unstable_matrix,
            p;
            dims=2,
        )

        plot_Δz_Δt(
            unstable_matrix,
            Δz,
            Δt,
            xlabel,
            ylabel,
            xscale=:log10,
            yscale=:log10,
            xticks=xticks,
            yticks=yticks,
            color=:black,
            a_I=a_I,
            fillalpha=0.5 * (1 - a_I),
        )
    end
    display(current())
end


function stability_check(p::SimulationParameters, n_zs, n_ts, var1, var2)
    domain = var1 == :n_A ? "atm" : "oce"
    unstable_matrix = zeros(length(n_ts), length(n_zs))

    for (i, n_z) in enumerate(n_zs)
        for (j, n_t) in enumerate(n_ts)
            setproperty!(p, var1, n_z)
            setproperty!(p, var2, n_t)

            cs = get_coupled_sim(p)
            sim = domain == "atm" ? cs.model_sims.atmos_sim : cs.model_sims.ocean_sim

            try
                Interfacer.step!(sim, p.Δt_cpl)
                unstable_matrix[i, j] = NaN
            catch err
                if isa(err, UnstableError)
                    unstable_matrix[i, j] = Inf
                else
                    rethrow()
                end
            end
        end
    end
    return unstable_matrix
end


function plot_ϱ_over_var(var_name; iterations=10, kwargs...)
    variable_dict = Dict(
        :Δt_cpl => [Base.logrange(10, 1e6, length=6), L"$\Delta t_{cpl}$"],
        :n_A => [Int.(Base.logrange(20, 2e5, length=5)), L"$\Delta z_A$"],
        :n_O => [Int.(Base.logrange(5, 5e5, length=6)), L"$\Delta z_O$"],
        :C_H_AI => [
            Float64.(LinRange(0.0008691059360985882, 0.0027815836784748733, 10)),
            L"$C_{H,AI}$",
        ],
        :C_H_AO => [
            Float64.(LinRange(0.0006765386884900067, 0.0014898726741724184, 10)),
            L"$C_{H,AO}$",
        ],
        :C_H_IO => [Float64.(LinRange(0.001, 0.01, 10)), L"$C_{H,OI}$"],
        :Δu_AO => [Float64.([0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]), L"$\Delta u$"],
        :n_t_A =>
            [Float64.(10 .^ LinRange(log10(1), log10(1000), 10)), L"$\Delta t_A$"],
        :n_t_O =>
            [Float64.(10 .^ LinRange(log10(1), log10(1000), 10)), L"$\Delta t_O$"],
        :a_I => [vec(0:0.1:1), L"a_I"],
    )
    var = variable_dict[var_name][1]
    par_name = variable_dict[var_name][2]
    physical_values = SimulationParameters(Δt_min=10, t_max=1000, Δt_cpl=1000, a_I=0.0, T_Ib=270.0)
    ϱs_atm, ϱs_oce, param_analytic, ϱs_analytic =
        get_ϱs_one_variable(
            physical_values,
            var,
            var_name,
            iterations=iterations,
            log_scale=false,
        )
    var, ϱs_oce, ϱs_atm, param_analytic, ϱs_analytic =
        treat_grid_sizes(
            var,
            var_name,
            ϱs_oce,
            ϱs_atm,
            physical_values;
            param_analytic=param_analytic,
            ϱs_analytic=ϱs_analytic,
        )
    gr()
    plot()
    # If there is instabilities, handle the text
    unstable_oce_indices = isinf.(ϱs_oce)
    ϱs_oce[unstable_oce_indices] .= NaN

    unstable_atm_indices = isinf.(ϱs_atm)
    ϱs_atm[unstable_atm_indices] .= NaN

    plot_data(
        ϱs_oce,
        ϱs_atm,
        var,
        par_name;
        ϱs_analytic=ϱs_analytic,
        param_analytic=param_analytic,
        color=:black,
        linestyles=:solid,
        yscale=:log10,
        kwargs...
    )
    display(current())
end


function get_ϱs_one_variable(
    p::SimulationParameters,
    vars,
    var_name;
    iterations=10,
    log_scale=false,
)
    ϱs_atm = zeros(length(vars))
    ϱs_oce = zeros(length(vars))

    for (k, var) in enumerate(vars)
        setproperty!(p, var_name, var)
        restore_physical_values!(p)
        if var_name == :Δt_cpl
            p.t_max = var
        end
        _, ϱs_atm[k], ϱs_oce[k] = run_simulation(p, iterations=iterations)
    end

    if log_scale
        finely_spaced_var = Base.logrange(vars[1], vars[end], length=100)
    else
        finely_spaced_var = range(vars[1], vars[end], length=100)
    end
    ϱs_analytic = zeros(length(finely_spaced_var))
    for (k, var) in enumerate(finely_spaced_var)
        setproperty!(p, var_name, var)
        restore_physical_values!(p)
        if var_name == :Δt_cpl
            p.t_max = var
        end
        ϱs_analytic[k] = compute_ϱ_analytical(p)
    end
    return ϱs_atm, ϱs_oce, finely_spaced_var, ϱs_analytic
end


function treat_grid_sizes(
    var,
    var_name,
    ϱs_oce,
    ϱs_atm,
    p::SimulationParameters;
    dims=1,
    param_analytic=nothing,
    ϱs_analytic=nothing,
)
    # Special treatment of n_A and n_O.
    if var_name == :n_A
        var = (p.h_A - p.z_A0) ./ reverse(var)
        if !isnothing(ϱs_analytic) && !isnothing(param_analytic)
            param_analytic = reverse(
                (p.h_A - p.z_A0) ./ param_analytic,
            )
        end
    elseif var_name == :n_O
        var = (p.h_O - p.z_O0) ./ reverse(var)
        if !isnothing(ϱs_analytic) && !isnothing(param_analytic)
            param_analytic = reverse(
                (p.h_O - p.z_O0) ./ param_analytic,
            )
        end
    elseif var_name == :n_t_A || var_name == :n_t_O
        var = p.Δt_min ./ reverse(var)
        if !isnothing(ϱs_analytic) && !isnothing(param_analytic)
            param_analytic = reverse(p.Δt_min ./ param_analytic)
        end
    end
    if var_name in [:n_A, :n_O, :n_t_A, :n_t_O]
        ϱs_oce =
            !isnothing(ϱs_oce) ? reverse(ϱs_oce, dims=dims) : nothing
        ϱs_atm =
            !isnothing(ϱs_atm) ? reverse(ϱs_atm, dims=dims) : nothing
        ϱs_analytic =
            !isnothing(ϱs_analytic) ? reverse(ϱs_analytic, dims=dims) :
            nothing
    end
    return var, ϱs_oce, ϱs_atm, param_analytic, ϱs_analytic
end
