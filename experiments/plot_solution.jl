using clima_playground
using Plots
using LaTeXStrings
import ClimaCore as CC

function load_timestep(cs; index=nothing)
    coords = []
    solution = []
    for model_sim in cs.model_sims
        u = model_sim.integrator.sol.u
        if isnothing(index)
            index = lastindex(u)
        end
        space = axes(u[index].data)
        coord_data = parent(CC.Spaces.coordinates_data(space).z)[:, 1]
        sol_data = vec(u[end].data)
        push!(coords, coord_data)
        push!(solution, sol_data)
    end
    return coords, solution
end

function plot_solution(; kwargs...)
    cs, _, _ = coupled_heat_equations(; kwargs...)
    coords, sols = load_timestep(cs)
    plot()
    colors = [:skyblue, :seagreen, :black]
    for (coord, sol, color) in zip(coords, sols, colors)
        plot!(sol, coord, color=color, linestyle=:solid, marker=:circle, legend=false, markerstrokecolor=color)
    end
    plot!(xlabel="Temperature [K]", ylabel="Height [m]", title=L"T(t_\mathrm{end})")
    display(current())
end
