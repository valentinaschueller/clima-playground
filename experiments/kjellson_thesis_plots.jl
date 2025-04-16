using clima_playground

function figure71()
    plot_ρ_over_a_i()
end

function figure72a()
    plot_ρ_over_var(10, "u_atm", a_is=[0.1, 0.4, 0.7])
end

function figure72b()
    plot_ρ_over_var(10, "u_oce", a_is=[0.1, 0.4, 0.7])
end

function figure73a()
    plot_obukhov_C_dependencies("C_AO")
end

function figure73b()
    plot_obukhov_C_dependencies("C_AI")
end

function figure74a()
    plot_ρ_over_var(10, "C_AO", a_is=[0.1, 0.4, 0.7])
end

function figure74b()
    plot_ρ_over_var(10, "C_AI", a_is=[0.1, 0.4, 0.7])
end

function figure74c()
    plot_ρ_over_var(10, "C_OI", a_is=[0.1, 0.4, 0.7])
end

function figure75()
    plot_ρ_over_var(10, "t_max", a_is=[0.1, 0.4, 0.7], xscale=:log10)
end

function figure76a()
    plot_ρ_over_var(10, "n_atm", a_is=[0.1, 0.4, 0.7], xscale=:log10)
end

function figure76b()
    plot_ρ_over_var(10, "n_oce", a_is=[0.1, 0.4, 0.7], xscale=:log10)
end

function figure77a()
    plot_unstable_range("atm", a_is=[0.1, 0.4, 0.7])
end

function figure77b()
    plot_unstable_range("oce", a_is=[0.1, 0.4, 0.7])
end

function ρ_dependence_ν_ω()
    analytical_convergence_factor_dependence()
end
