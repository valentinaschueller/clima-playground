module clima_playground

export coupled_heat_equations
export plot_obukhov_C_dependencies
export analytical_convergence_factor_dependence
export compute_ρ_analytical
export define_realistic_vals

using LinearAlgebra
using Plots
using Statistics
using LaTeXStrings

include("heat_equations.jl")
include("analytical_plots.jl")

end
