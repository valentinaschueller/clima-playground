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

function plot_ice_thickness_convergence(; iterations=10, kwargs...)
    p = SimulationParameters(t_max=1000, Δt_cpl=1000, a_I=1.0, ice_model_type=:temp_feedback, C_H_IO=7.5e-3)
    h_Is = Base.logrange(1e-4, 1e3, length=7)
    ϱs_atm = zeros(length(h_Is))
    for (k, h_I) in enumerate(h_Is)
        setproperty!(p, :h_I_ini, h_I)
        _, ϱ_atm, ϱ_oce = run_simulation(p, iterations=iterations)
        ϱs_atm[k], _ = extract_ρ(ϱ_atm, ϱ_oce)
    end
    gr()
    plot()

    unstable_atm_indices = isinf.(ϱs_atm)
    ϱs_atm[unstable_atm_indices] .= NaN

    plot!(
        h_Is,
        ϱs_atm;
        label=L"$ϱ_\mathrm{num}$",
        markershape=:x,
        color=:black,
        linewidth=2,
        legendfontsize=12,
        xlabel=L"h_I",
        xscale=:log10,
        yscale=:log10,
        ylabel="ϱ",
        legend=:right,
        kwargs...
    )
    display(current())
end
