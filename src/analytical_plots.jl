using Plots


function plot_C_AO_dependence()
    params = SimulationParameters()

    neg_vars = -200:-50
    pos_vars = 10:200

    C_neg = zeros(length(neg_vars))
    C_pos = zeros(length(pos_vars))
    for (j, neg_var) in enumerate(neg_vars)
        C_neg[j] = compute_C_H_AO(params; L_AO=neg_var)
    end
    for (j, pos_var) in enumerate(pos_vars)
        C_pos[j] = compute_C_H_AO(params; L_AO=pos_var)
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
    params = SimulationParameters()
    a_is = 0:0.01:1

    neg_vars = -200:-50
    pos_vars = 10:200
    C_neg = zeros(length(a_is), length(neg_vars))
    C_pos = zeros(length(a_is), length(pos_vars))
    for (i, a_i) in enumerate(a_is)
        params.a_i = a_i
        for (j, neg_var) in enumerate(neg_vars)
            C_neg[i, j] = compute_C_H_AI(params; L_AI=neg_var)
        end
        for (j, pos_var) in enumerate(pos_vars)
            C_pos[i, j] = compute_C_H_AI(params; L_AI=pos_var)
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
    params = SimulationParameters()

    # Loop over values
    for (k, a_i) in enumerate(a_is)
        for (i, ν) in enumerate(νs)
            for (j, ω) in enumerate(ωs)
                params.a_i = a_i
                ρ = compute_ρ_analytical_with_ν_and_ω_variable(params, ν, ω)
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


function compute_ρ_analytical_with_ν_and_ω_variable(params, ν, ω)
    real_part_o =
        sqrt((sqrt(ν^2 + ω^2) + ν) / (2 * params.α_o))
    imag_part_o =
        im *
        sign(ω) *
        sqrt((sqrt(ν^2 + ω^2) - ν) / (2 * params.α_o))
    σ_o = real_part_o + imag_part_o
    real_part_a =
        sqrt((sqrt(ν^2 + ω^2) + ν) / (2 * params.α_a))
    imag_part_a =
        im *
        sign(ω) *
        sqrt((sqrt(ν^2 + ω^2) - ν) / (2 * params.α_a))
    σ_a = real_part_a + imag_part_a
    return abs(
        (1 - params.a_i)^2 * params.C_AO^2 / (
            (
                params.k_oce *
                σ_o *
                (1 / tanh(σ_o * (params.h_oce - params.z_0numO))) +
                (1 - params.a_i) * params.C_AO +
                params.a_i * params.C_IO
            ) * (
                params.k_atm *
                σ_a *
                (1 / tanh(σ_a * (params.h_atm - params.z_0numA))) +
                (1 - params.a_i) * params.C_AO +
                params.a_i * params.C_AI
            )
        ),
    )
end
