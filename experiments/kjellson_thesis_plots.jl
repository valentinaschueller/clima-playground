using clima_playground

function figure71()
    coupled_heat_equations(iterations=10, a_is=0:0.1:1, yscale=:log10, analytic_conv_fac=true, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000))
end

function figure72a()
    coupled_heat_equations(iterations=10, var_name="u_atm", a_is=[0.1, 0.4, 0.7], analytic_conv_fac=true, yscale=:log10, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000), compute_oce_conv_fac=false)
end

function figure72b()
    coupled_heat_equations(iterations=10, var_name="u_oce", a_is=[0.1, 0.4, 0.7], analytic_conv_fac=true, yscale=:log10, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000), compute_oce_conv_fac=false)
end

function figure73a()
    plot_obukhov_C_dependencies("C_AO")
end

function figure73b()
    plot_obukhov_C_dependencies("C_AI")
end

function figure74a()
    coupled_heat_equations(iterations=10, var_name="C_AO", a_is=[0.1, 0.4, 0.7], analytic_conv_fac=true, yscale=:log10, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000), compute_oce_conv_fac=false)
end

function figure74b()
    coupled_heat_equations(iterations=10, var_name="C_AI", a_is=[0.1, 0.4, 0.7], analytic_conv_fac=true, yscale=:log10, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000), compute_oce_conv_fac=false)
end

function figure74c()
    coupled_heat_equations(iterations=10, var_name="C_OI", a_is=[0.1, 0.4, 0.7], analytic_conv_fac=true, yscale=:log10, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000), compute_oce_conv_fac=false)
end

function figure75()
    coupled_heat_equations(iterations=10, var_name="t_max", a_is=[0.1, 0.4, 0.7], analytic_conv_fac=true, yscale=:log10, xscale=:log10, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000), compute_oce_conv_fac=false)
end

function figure76a()
    coupled_heat_equations(iterations=10, var_name="n_atm", a_is=[0.1, 0.4, 0.7], analytic_conv_fac=true, yscale=:log10, xscale=:log10, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000), compute_oce_conv_fac=false)
end

function figure76b()
    coupled_heat_equations(iterations=10, var_name="n_oce", a_is=[0.1, 0.4, 0.7], analytic_conv_fac=true, yscale=:log10, xscale=:log10, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000), compute_oce_conv_fac=false)
end

function figure77a()
    plot_unstable_range("atm", a_is=[0.1, 0.4, 0.7])
end

function figure77b()
    plot_unstable_range("oce", a_is=[0.1, 0.4, 0.7])
end
