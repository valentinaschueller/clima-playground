module clima_playground

export coupled_heat_equations, solve_coupler!, run_simulation

include("parameters.jl")
include("analysis.jl")
include("diagnostics.jl")
include("components/atmosphere.jl")
include("components/ocean.jl")
include("components/ice.jl")
include("coupled_simulation.jl")
include("monin_obukhov.jl")
include("postprocessing.jl")

function run_simulation(
    physical_values;
    iterations=10,
    parallel=false,
)
    cs = get_coupled_sim(physical_values)
    ϱ_A, ϱ_O = solve_coupler!(
        cs,
        iterations=iterations,
        parallel=parallel,
    )
    return cs, ϱ_A, ϱ_O
end

function coupled_heat_equations(;
    iterations::Int=1,
    parallel::Bool=false,
    monin_obukhov::Bool=true,
    kwargs...,
)
    physical_values = SimulationParameters{Float64}(; kwargs...)

    if monin_obukhov
        physical_values.C_H_AO = compute_C_H_AO(physical_values)
        physical_values.C_H_AI = compute_C_H_AI(physical_values)
        restore_physical_values!(physical_values)
    end

    cs, ϱ_A, ϱ_O = run_simulation(physical_values, iterations=iterations, parallel=parallel)
    return cs, ϱ_A, ϱ_O
end;


end
