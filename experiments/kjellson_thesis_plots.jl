using clima_playground
using Plots
using LaTeXStrings

function figure71()
    plot_ρ_over_a_i()
end

function figure72a()
    plot_ρ_over_var(10, :Δu_AO, a_is=[0.1, 0.4, 0.7])
end

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

function figure74a()
    plot_ρ_over_var(10, :C_H_AO, a_is=[0.1, 0.4, 0.7])
end

function figure74b()
    plot_ρ_over_var(10, :C_H_AI, a_is=[0.1, 0.4, 0.7])
end

function figure74c()
    plot_ρ_over_var(10, :C_H_IO, a_is=[0.1, 0.4, 0.7])
end

function figure75()
    plot_ρ_over_var(10, :Δt_cpl, a_is=[0.1, 0.4, 0.7], xscale=:log10)
end

function figure76a()
    plot_ρ_over_var(10, :n_atm, a_is=[0.1, 0.4, 0.7], xscale=:log10)
end

function figure76b()
    plot_ρ_over_var(10, :n_oce, a_is=[0.1, 0.4, 0.7], xscale=:log10)
end

function figure77()
    plot_unstable_range("atm", a_is=[0.1, 0.4, 0.7])
end

function figure78()
    plot_unstable_range("oce", a_is=[0.1, 0.4, 0.7])
end

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
                ρ = compute_ρ_analytical(params; s=ν + im * ω)
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
