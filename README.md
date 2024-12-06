# Toy Examples for ClimaCoupler.jl

⚠️ This is work in progress and not running at the moment! ⚠️

This repository contains different coupled (toy) problems and coupling setups that are representative of and relevant for atmosphere-ocean-sea ice coupling.
The idea is to explicitly make use of the [ClimaCoupler.jl](https://github.com/CliMA/ClimaCoupler.jl) API for these problems.

To run the code:
```bash
> julia --project
> ]instantiate
# return to REPL
> using clima_playground
> coupled_heat_equations()
```
