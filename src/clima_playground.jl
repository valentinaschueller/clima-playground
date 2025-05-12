module clima_playground

export coupled_heat_equations
export run_simulation
export compute_ϱ_analytical
export SimulationParameters
export compute_C_H_AI
export compute_C_H_AO
export restore_physical_values!
export extract_ρ
export stability_check
export solve_surface_energy_balance

using LinearAlgebra
using Plots
using Statistics

include("heat_equations.jl")

end
