using clima_playground
using Plots
using LaTeXStrings
using NetCDF
import ClimaCore as CC
import ClimaDiagnostics as CD

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
    p = cs.model_sims[1].params
    plot()
    colors = [:skyblue, :seagreen, :black]
    if p.a_I == 0.0
        coords, sols, colors = coords[1:2], sols[1:2], colors[1:2]
    end
    for (coord, sol, color) in zip(coords, sols, colors)
        plot!(sol, coord, color=color, linestyle=:solid, marker=:circle, legend=false, markerstrokecolor=color)
    end
    plot!(xlabel="Temperature [K]", ylabel="Height [m]", title=L"T(t_\mathrm{end})")
    display(current())
end

function plot_solution_over_time(; kwargs...)
    cs, _, _ = coupled_heat_equations(; kwargs...)
    p = cs.model_sims[1].params
    dt = CD.seconds_to_str_short(p.Δt_min)
    time = ncread("output/h_I_$(dt)_inst.nc", "time")
    T_O = ncread("output/T_O_$(dt)_inst.nc", "T_O", start=[1, p.n_O], count=[-1, 1])
    T_A = ncread("output/T_A_$(dt)_inst.nc", "T_A", start=[1, 1], count=[-1, 1])
    T_Is = ncread("output/T_Is_$(dt)_inst.nc", "T_Is")
    h_I = ncread("output/h_I_$(dt)_inst.nc", "h_I")
    p1 = plot(time, [T_A T_O T_Is], xlabel="Time [s]", ylabel="Temperature [K]", label=[L"T_A" L"T_O" L"T_{I,s}"], color=[:skyblue :seagreen :black])
    p2 = plot(time, h_I, color=:black, label=L"h_I", ylabel="Ice Thickness [m]", xlabel="Time [s]")
    l = @layout [a b]
    plot(p1, p2, layout=l, legendfontsize=12, linewidth=2)
end
