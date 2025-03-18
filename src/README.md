To get the same behaviour as in the main branch, run:
coupled_heat_equations(iterate=1, values=Dict(:delta_t_min=>1, :t_max=>3600, :delta_t_cpl=>100))

The times are specified explicitly here as default is delta_t_min=10 and t_max=delta_t_cpl=1000, 
with equal t_max and delta_t_cpl to be able to only compute the convergence factor once.

To get a plot of convergence factor as a function of iteration, run:
coupled_heat_equations(plot_conv_facs_iter=true)

If you want to use the parallel Schwarz iteration, run:
coupled_heat_equations(plot_conv_facs_iter=true, parallel=true)

To change to closest in time boundary mapping, run:
coupled_heat_equations(plot_conv_facs_iter=true, parallel=true, boundary_mapping="cit")

To combine two on eachother following convergence factors to see that the parallel and alternating Schwarz iteration are related, run:
coupled_heat_equations(print_conv_facs_iter=true, atm=true, parallel=true, boundary_mapping="cit", combine=true)

and for the alternating Schwarz:
coupled_heat_equations(print_conv_facs_iter=true, atm=true, parallel=false, boundary_mapping="cit")

To plot the convergence factor wrt a_i, run:
coupled_heat_equations(a_is=0:0.1:1, log_conv_fac=true, analytic_conv_fac=true)

To plot the convergence factor as a function of another variable, run:
coupled_heat_equations(var_name="u_atm", values=Dict(:a_i=>0.1), analytic_conv_fac=true)

For the same plot but for different values of a_i, run:
coupled_heat_equations(var_name="u_atm", a_is=[0.1, 0.4, 0.7], analytic_conv_fac=true)

To get a plot of the cfl condition, run:
coupled_heat_equations(plot_unstable_range=true, a_is=[0.1])
This also plots a theoretical limit for the cfl condition, based on the maximum of the terms in the boundary update.

To plot the same limit for different a_i, run:
coupled_heat_equations(plot_unstable_range=true, a_is=[0.1, 0.4, 0.7])




coupled_heat_equations(; iterate=10, parallel=false, boundary_mapping="mean", values=Dict{Symbol,Int}(), print_conv_facs_iter=false, plot_conv_facs_iter=false, analytic_conv_fac=false, atm=false, combine=false, plot_unstable_range=false, a_is=[], var_name=nothing, xscale = :identity, legend = :topright, log_conv_fac=false, xticks=nothing)