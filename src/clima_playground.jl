module clima_playground

export coupled_heat_equations
export plot_C_AI_dependence
export plot_C_AO_dependence
export analytical_convergence_factor_dependence
export compute_ρ_analytical
export define_realistic_vals
export plot_ρ_over_k
export plot_unstable_range
export plot_ρ_over_a_i
export plot_ρ_over_var

using LinearAlgebra
using Plots
using Statistics
using LaTeXStrings

include("heat_equations.jl")
include("analytical_plots.jl")

end
