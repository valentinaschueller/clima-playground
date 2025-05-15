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
        sol_data = vec(u[index].data)
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

function plot_solution_over_time(; kwargs...)
    cs, _, _ = coupled_heat_equations(; kwargs...)
    time = cs.model_sims.atmos_sim.integrator.sol.t
    T_A = get_field(cs.model_sims.atmos_sim, Val(:T_atm_sfc))
    T_O = get_field(cs.model_sims.ocean_sim, Val(:T_oce_sfc))
    T_I = get_field(cs.model_sims.ice_sim, Val(:T_ice))
    h_I = get_field(cs.model_sims.ice_sim, Val(:h_I))
    p1 = plot(time, [T_A T_O T_I], xlabel="Time [s]", ylabel="Temperature [K]", label=[L"T_A" L"T_O" L"T_I"], color=[:skyblue :seagreen :black])
    p2 = plot(time, h_I, color=:black, label=L"h_I", ylabel="Ice Thickness [m]", xlabel="Time [s]")
    l = @layout [a b]
    plot(p1, p2, layout=l, legendfontsize=12, linewidth=2)
end
