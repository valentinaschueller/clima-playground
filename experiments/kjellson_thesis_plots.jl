using clima_playground

function figure71()
    plot_ρ_over_a_i()
end

function figure72a()
    plot_ρ_over_var(10, :Δu_AO, a_is=[0.1, 0.4, 0.7])
end

function figure73a()
    plot_C_AO_dependence()
end

function figure73b()
    plot_C_AI_dependence()
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

function ρ_dependence_ν_ω()
    analytical_convergence_factor_dependence()
end
