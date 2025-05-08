module clima_playground

export coupled_heat_equations
export compute_ρ_analytical
export SimulationParameters
export compute_C_H_AI
export compute_C_H_AO
export restore_physical_values!
export plot_ρ_over_k
export plot_unstable_range
export plot_ρ_over_a_i
export plot_ρ_over_var

using LinearAlgebra
using Plots
using Statistics
using LaTeXStrings

include("heat_equations.jl")

end
