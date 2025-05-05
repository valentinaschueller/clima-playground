using Plots


function plot_C_AO_dependence()
    physical_values = define_realistic_vals()

    neg_vars = -200:-50
    pos_vars = 10:200

    C_neg = zeros(length(neg_vars))
    C_pos = zeros(length(pos_vars))
    for (j, neg_var) in enumerate(neg_vars)
        physical_values[:L_AO] = neg_var
        physical_values[:L_AI] = physical_values[:L_AI]
        compute_C_AO!(physical_values)
        C_neg[j] = physical_values[:C_AO]
    end
    for (j, pos_var) in enumerate(pos_vars)
        physical_values[:L_AO] = pos_var
        physical_values[:L_AI] = physical_values[:L_AI]
        compute_C_AO!(physical_values)
        C_pos[j] = physical_values[:C_AO]
    end

    println(maximum([maximum(C_neg), maximum(C_pos)]))
    println(minimum([minimum(C_neg), minimum(C_pos)]))

    gr()
    plot(
        neg_vars,
        C_neg,
        xlabel=L"$L^A_O$",
        ylabel=L"$C^A_O$",
        label="",
        color=:black,
    )
    plot!(pos_vars, C_pos, label="", color=:black)
    display(current())
end


function plot_C_AI_dependence()
    physical_values = define_realistic_vals()
    a_is = 0:0.01:1

    neg_vars = -200:-50
    pos_vars = 10:200
    C_neg = zeros(length(a_is), length(neg_vars))
    C_pos = zeros(length(a_is), length(pos_vars))
    for (i, a_i) in enumerate(a_is)
        physical_values[:a_i] = a_i
        for (j, neg_var) in enumerate(neg_vars)
            physical_values[:L_AI] = neg_var
            C_neg[i, j] = compute_C_AI(physical_values)
        end
        for (j, pos_var) in enumerate(pos_vars)
            physical_values[:L_AI] = pos_var
            C_pos[i, j] = compute_C_AI(physical_values)
        end
    end

    println(maximum([maximum(C_neg), maximum(C_pos)]))
    println(minimum([minimum(C_neg), minimum(C_pos)]))

    plotly()
    surface(a_is, neg_vars, C_neg', xlabel="aᴵ", ylabel="Lᴬᴵ", zlabel="Cᴬᴵ")
    surface!(a_is, pos_vars, C_pos')
end


"""Computes and plots the analytical convergence factor as a function of `ν` and `ω`."""
function analytical_convergence_factor_dependence()
    νs = range(0, stop=10, length=50)
    ωs = range(0.001, stop=10, length=50)
    a_is = range(0.001, stop=1, length=10)
    ρs = zeros(length(νs), length(ωs), length(a_is))
    physical_values = define_realistic_vals()

    # Loop over values
    for (k, a_i) in enumerate(a_is)
        for (i, ν) in enumerate(νs)
            for (j, ω) in enumerate(ωs)
                physical_values[:ω] = ω
                physical_values[:ν] = ν
                physical_values[:a_i] = a_i
                ρ = compute_ρ_analytical_with_ν_and_ω_variable(physical_values)
                ρs[i, j, k] = ρ
            end
        end
        i_max, j_max = Tuple(CartesianIndices(ρs)[argmax(ρs)])
        ν_max = νs[i_max]
        ω_max = ωs[j_max]
        println("supremum at ν=$ν_max and ω=$ω_max")
    end


    plotly()
    surface(
        ωs,
        νs,
        ρs[:, :, 1],
        color=:viridis,
        xlabel="ω",
        ylabel="ν",
        zlabel="̂ρ(ν+iω)",
    )
end

"""
Computes the analytical convergence factor as a function of `ν` and `ω`.

**Arguments:**

-`params::Dict`: Can be defined using `define_realistic_vals()`, but should also have the keys `:omega` and `:nu`.

"""
function compute_ρ_analytical_with_ν_and_ω_variable(params)
    real_part_o =
        sqrt((sqrt(params[:ν]^2 + params[:ω]^2) + params[:ν]) / (2 * params[:α_o]))
    imag_part_o =
        im *
        sign(params[:ω]) *
        sqrt((sqrt(params[:ν]^2 + params[:ω]^2) - params[:ν]) / (2 * params[:α_o]))
    σ_o = real_part_o + imag_part_o
    real_part_a =
        sqrt((sqrt(params[:ν]^2 + params[:ω]^2) + params[:ν]) / (2 * params[:α_a]))
    imag_part_a =
        im *
        sign(params[:ω]) *
        sqrt((sqrt(params[:ν]^2 + params[:ω]^2) - params[:ν]) / (2 * params[:α_a]))
    σ_a = real_part_a + imag_part_a
    η_AO =
        params[:C_AO] *
        abs(params[:u_atm] - params[:u_oce]) *
        params[:ρ_atm] *
        params[:c_atm]
    η_OI = params[:C_OI] * abs(params[:u_oce]) * params[:ρ_oce] * params[:c_oce]
    η_AI = params[:C_AI] * abs(params[:u_atm]) * params[:ρ_atm] * params[:c_atm]
    return abs(
        (1 - params[:a_i])^2 * η_AO^2 / (
            (
                params[:k_oce] *
                σ_o *
                (1 / tanh(σ_o * (params[:h_oce] - params[:z_0numO]))) +
                (1 - params[:a_i]) * η_AO +
                params[:a_i] * η_OI
            ) * (
                params[:k_atm] *
                σ_a *
                (1 / tanh(σ_a * (params[:h_atm] - params[:z_0numA]))) +
                (1 - params[:a_i]) * η_AO +
                params[:a_i] * η_AI
            )
        ),
    )
end
