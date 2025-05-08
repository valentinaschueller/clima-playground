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

Simulation parameters (physical and numerical) can be changed by passing keyword arguments, e.g.,
```julia
coupled_heat_equations(iterations=10, Δt_min=10, t_max=1000, Δt_cpl=1000, boundary_mapping="mean")
```

If you want to use the parallel Schwarz iteration, run:
```julia
coupled_heat_equations(iterations=10, parallel=true)
```

To reproduce the figures from Hanna Kjellson's MSc thesis, run, e.g.:
```julia
include("experiments/kjellson_thesis_plots.jl")
figure71() # creates the plot corresponding to Fig. 7.1 in the thesis
```
