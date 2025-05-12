using clima_playground
using Plots
using LaTeXStrings


function plot_ice_seb_results()
    params = SimulationParameters(C_H_AI=1.4e-3)
    T_As = vec(range(260, 290, length=100))
    caches = [(T_A=T_A, C_AI=params.C_AI) for T_A in T_As]
    T_ice = solve_surface_energy_balance.(caches)
    plot(T_As, T_ice, color=:black, xlabel=L"T_A", ylabel=L"T_{I,s}", legend=false)
    display(current())
end
