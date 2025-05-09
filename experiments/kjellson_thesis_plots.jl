using clima_playground
using Plots
using LaTeXStrings

function figure71()
    plot_ρ_over_a_i()
end

function figure72()
    plot_ρ_over_var(10, :Δu_AO, a_is=[0.1, 0.4, 0.7])
end

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
    a_is = 0:0.01:1
    L_AIs = vec(-200:200)
    params = SimulationParameters()
    C_H_AI = zeros(length(a_is), length(L_AIs))
    for (i, a_i) in enumerate(a_is)
        params.a_i = a_i
        for (j, L_AI) in enumerate(L_AIs)
            C_H_AI[i, j] = compute_C_H_AI(params; L_AI=L_AI)
        end
    end
    plot()
    surface(L_AIs, a_is, C_H_AI, xlabel=L"L_{AI}", ylabel=L"a_I", zlabel=L"C_{H,AI}")
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
