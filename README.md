# Toy Examples for ClimaCoupler.jl

⚠️ This is work in progress! ⚠️

This repository contains coupled (toy) problems and coupling setups that are representative of and relevant for atmosphere-ocean-sea ice coupling.
The idea is to explicitly make use of the [ClimaCoupler.jl](https://github.com/CliMA/ClimaCoupler.jl) API for these problems.

To run the code:
```bash
> julia --project
> ]instantiate
# return to REPL
> using clima_playground
> coupled_heat_equations()
```
This solves the coupled heat equation on the time interval $t\in[0:3600] s$ with coupling time step $\Delta t_{cpl}=100 s$ and only one iteration per coupling time step.

To get a plot of the convergence factor as a function of iteration, run:
```julia
coupled_heat_equations(iterations=10, plot_conv_facs_iter=true, params=Dict(:Δt_min=>10, :t_max=>1000, :Δt_cpl=>1000))
```

If you want to use the parallel Schwarz iteration, run:
```julia
coupled_heat_equations(iterations=10, plot_conv_facs_iter=true, parallel=true, params=Dict(:Δt_min=>10, :t_max=>1000, :Δt_cpl=>1000))
```

To change to closest in time boundary mapping, run:
```julia
coupled_heat_equations(iterations=10, plot_conv_facs_iter=true, parallel=true, boundary_mapping="cit", params=Dict(:Δt_min=>10, :t_max=>1000, :Δt_cpl=>1000))
```

To combine two on eachother following convergence factors to see that the parallel and alternating Schwarz iteration are related, run:
```julia
coupled_heat_equations(iterations=10, plot_conv_facs_iter=true, parallel=true, boundary_mapping="cit", combine_ρ_parallel=true, params=Dict(:Δt_min=>10, :t_max=>1000, :Δt_cpl=>1000))
```

and for the alternating Schwarz:
```julia
coupled_heat_equations(iterations=10, plot_conv_facs_iter=true, parallel=false, boundary_mapping="cit", params=Dict(:Δt_min=>10, :t_max=>1000, :Δt_cpl=>1000))
```

To plot the convergence factor wrt a_i, run:
```julia
coupled_heat_equations(iterations=10, a_is=0:0.1:1, yscale=:log10, analytic_conv_fac=true, params=Dict(:Δt_min=>10, :t_max=>1000, :Δt_cpl=>1000))
```

To plot the convergence factor as a function of another variable, run for instance:
```julia
coupled_heat_equations(iterations=10, var_name="u_atm", analytic_conv_fac=true, yscale=:log10, params=Dict(:a_i=>0.1, :Δt_min=>10, :t_max=>1000, :Δt_cpl=>1000))
```

For the same plot but for different values of a_i, run:
```julia
coupled_heat_equations(iterations=10, var_name="u_atm", a_is=[0.1, 0.4, 0.7], analytic_conv_fac=true, yscale=:log10, params=Dict(:Δt_min=>10, :t_max=>1000, :Δt_cpl=>1000))
```

To get a plot of the cfl condition for the atmosphere, run:
```julia
coupled_heat_equations(plot_unstable_range=true, a_is=[0.1], params=Dict(:Δt_min=>10, :t_max=>1000, :Δt_cpl=>1000))
```

To plot the same limit for three different a_i, run:
```julia
coupled_heat_equations(plot_unstable_range=true, a_is=[0.1, 0.4, 0.7], params=Dict(:Δt_min=>10, :t_max=>1000, :Δt_cpl=>1000))
```

To get the corresponding plot for the ocean, run:
```julia
coupled_heat_equations(plot_unstable_range=true, a_is=[0.1, 0.4, 0.7], params=Dict(:Δt_min=>10, :t_max=>1000, :Δt_cpl=>1000), compute_atm_conv_fac=false)
```

To verify that the convergence factor is a decreasing function of $\omega$ and $\nu$, run:
```julia
analytical_convergence_factor_dependence()
```

 To plot the dependence of C_AO on L_AO, run:
 ```julia 
 plot_obukhov_C_dependencies("C_AO")
 ```

 To plot the dependence of C_AI on L_AI and a_i, run:
 ```julia
 plot_obukhov_C_dependencies("C_AI")
 ```
