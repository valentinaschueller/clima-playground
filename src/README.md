To solve the coupled heat equation on the time interval $t\in[0:3600] s$ with coupling time step $\Delta t_{cpl}=100 s$ and only one iteration per coupling time step, run:
```julia
coupled_heat_equations()
```

To get a plot of the convergence factor as a function of iteration, run:
```julia
coupled_heat_equations(iterate=10, plot_conv_facs_iter=true, values=Dict(:delta_t_min=>10, :t_max=>1000, :delta_t_cpl=>1000))
```

If you want to use the parallel Schwarz iteration, run:
```julia
coupled_heat_equations(iterate=10, plot_conv_facs_iter=true, parallel=true, values=Dict(:delta_t_min=>10, :t_max=>1000, :delta_t_cpl=>1000))
```

To change to closest in time boundary mapping, run:
```julia
coupled_heat_equations(iterate=10, plot_conv_facs_iter=true, parallel=true, boundary_mapping="cit", values=Dict(:delta_t_min=>10, :t_max=>1000, :delta_t_cpl=>1000))
```

To combine two on eachother following convergence factors to see that the parallel and alternating Schwarz iteration are related, run:
```julia
coupled_heat_equations(iterate=10, print_conv_facs_iter=true, parallel=true, boundary_mapping="cit", combine=true, values=Dict(:delta_t_min=>10, :t_max=>1000, :delta_t_cpl=>1000))
```

and for the alternating Schwarz:
```julia
coupled_heat_equations(iterate=10, print_conv_facs_iter=true, parallel=false, boundary_mapping="cit", values=Dict(:delta_t_min=>10, :t_max=>1000, :delta_t_cpl=>1000))
```

To plot the convergence factor wrt a_i, run:
```julia
coupled_heat_equations(iterate=10, a_is=0:0.1:1, yscale=:log10, analytic_conv_fac=true, values=Dict(:delta_t_min=>10, :t_max=>1000, :delta_t_cpl=>1000))
```

To plot the convergence factor as a function of another variable, run:
```julia
coupled_heat_equations(iterate=10, var_name="u_atm", values=Dict(:a_i=>0.1), analytic_conv_fac=true, yscale=:log10, values=Dict(:delta_t_min=>10, :t_max=>1000, :delta_t_cpl=>1000))
```

For the same plot but for different values of a_i, run:
```julia
coupled_heat_equations(iterate=10, var_name="u_atm", a_is=[0.1, 0.4, 0.7], analytic_conv_fac=true, yscale=:log10, values=Dict(:delta_t_min=>10, :t_max=>1000, :delta_t_cpl=>1000))
```

To get a plot of the cfl condition, run:
```julia
coupled_heat_equations(plot_unstable_range=true, a_is=[0.1], values=Dict(:delta_t_min=>10, :t_max=>1000, :delta_t_cpl=>1000))
```

To plot the same limit for three different a_i, run:
```julia
coupled_heat_equations(plot_unstable_range=true, a_is=[0.1, 0.4, 0.7], values=Dict(:delta_t_min=>10, :t_max=>1000, :delta_t_cpl=>1000))
``