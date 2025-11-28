# Toy Examples for ClimaCoupler.jl

---
⚠️ This branch is to experiment with using the `ClimaAtmos.jl` SCM as an atmosphere component ⚠️

You need to locally check out the latest versions of `ClimaCoupler.jl` (I am at [`a3dc89e`](https://github.com/CliMA/ClimaCoupler.jl/commit/a3dc89e)), as well as the correct branches for the atmospheric SCM:
- `ClimaAtmos.jl`: [`jy/larcform1_1M`](https://github.com/CliMA/ClimaAtmos.jl/tree/jy/larcform1_1M)
- `AtmosphericProfilesLibrary.jl`: [`jy/Larcform1`](https://github.com/CliMA/AtmosphericProfilesLibrary.jl/tree/jy/Larcform1)

You should link to them when setting up `clima_playground` (see below)

```julia
> ]dev /path/to/ClimaCoupler.jl /path/to/ClimaAtmos.jl /path/to/AtmosphericProfilesLibrary.jl
> instantiate
> using clima_playground
> coupled_run(config_file="clima-playground/experiments/config.yaml", reference_job_id="single_column_precipitation_test")
```

Relevant files to look at:
- `src/coupled_driver.jl`: contains `atmos_only` and `coupled_run`, to run a standalone atmosphere or coupled simulation, respectively. Runs the model, then creates reference plots as for `ClimaAtmos.jl` runs.
  - `atmos_only(config_file="path/to/config/file.yaml", reference_job_id="ref_job_id")`: should create the same results as a ClimaAtmos run with the same config file
  - `coupled_run(config_file="path/to/config/file.yaml", reference_job_id="ref_job_id")`: runs the atmosphere, coupled to a surface with constant surface temperature equal to the value provided in the config file.
- `src/components/climaatmos_scm.jl`: the atmosphere component. Essentially a copy of [`climaatmos.jl`](https://github.com/CliMA/ClimaCoupler.jl/blob/main/experiments/ClimaEarth/components/atmosphere/climaatmos.jl) in the `ClimaEarth` examples. Crucial difference: I had to add a specific instantiation of `Interfacer.remap()` since remapping is a lot simpler for the single column setup
- `src/components/simple_surface.jl`: an "ice model" that really only returns a constant surface temperature, provided in the config file with key `T_sfc`. I have put everything in place s.t. one could introduce a dynamic ice thickness and make the surface temperature model more advanced.
  - right now the incoming atmospheric fields are not actually used but this would be the way to go forward.
- as a config file I have provided `experiments/config.yaml`, which is a copy of [`single_column_precipitation_test.yaml`](https://github.com/CliMA/ClimaAtmos.jl/blob/main/config/model_configs/single_column_precipitation_test.yml) with some extra fields used by the surface model (`T_sfc`, `alb_I`,...). Remember to check the path to the TOML!

---

This repository contains coupled (toy) problems and coupling setups that are representative of and relevant for atmosphere-ocean-ice coupling.
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
coupled_heat_equations(iterations=10, Δt_min=10, t_max=1000, Δt_cpl=1000)
```

If you want to use the parallel Schwarz iteration, run:
```julia
coupled_heat_equations(iterations=10, parallel=true)
```

To activate sea ice, use a nonzero ice area fraction, e.g., `a_I=0.5`.
There are three implemented ice models:
- `:constant`: ice thickness and temperature remain constant
- `:temperature_feedback`: constant ice thickness, dynamic ice temperature
- `:thickness_feedback`: dynamic ice thickness and temperature (the 0-layer model from Semtner, [1976](https://doi.org/10.1175/1520-0485(1976)006<0379:AMFTTG>2.0.CO;2))

The directory `experiments` contains scripts for plotting analytical and numerical results with this code.
- `kjellson_thesis_plots.jl` creates plots resembling figures from Hanna Kjellson's MSc thesis. 
- `plot_solution.jl` creates plots that show the solution of a simulation over time or space
- `plot_ice.jl` contains plots for when the dynamic sea ice component is active
- `semtner_testcase.jl` runs a validation test case for the ice model based on (Semtner, [1976](https://doi.org/10.1175/1520-0485(1976)006<0379:AMFTTG>2.0.CO;2))
