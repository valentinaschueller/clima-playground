using clima_playground
using Plots
using LaTeXStrings
import ClimaCore as CC

function plot_ice_seb_results()
    params = SimulationParameters(C_H_AI=1.4e-3)
    T_As = vec(range(260, 290, length=100))
    caches = [(T_A=T_A, C_AI=params.C_AI) for T_A in T_As]
    T_ice = solve_surface_energy_balance.(caches)
    plot(T_As, T_ice, color=:black, xlabel=L"T_A", ylabel=L"T_{I,s}", legend=false)
    display(current())
end

function plot_ice_thickness_convergence(; iterations=5, kwargs...)
    p = SimulationParameters(Δt_min=200, t_max=3600, Δt_cpl=3600, a_I=1.0, ice_model_type=:thickness_feedback)
    h_Is = Base.logrange(1e-4, 1e3, length=10)
    ϱs_atm = zeros(length(h_Is))
    for (k, h_I) in enumerate(h_Is)
        setproperty!(p, :h_I_ini, h_I)
        _, ϱs_atm[k], _ = run_simulation(p, iterations=iterations)
    end

    plot(
        h_Is,
        ϱs_atm;
        label=L"$ϱ_\mathrm{num}, h=h(t)$",
        markershape=:o,
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

    p.ice_model_type = :temp_feedback
    ϱs_atm = zeros(length(h_Is))
    for (k, h_I) in enumerate(h_Is)
        setproperty!(p, :h_I_ini, h_I)
        _, ϱs_atm[k], _ = run_simulation(p, iterations=iterations)
    end

    plot!(
        h_Is,
        ϱs_atm;
        label=L"$ϱ_\mathrm{num}, h=h_0$",
        markershape=:x,
        color=:black,
        linestyle=:dash,
        linewidth=2,
        legendfontsize=12,
        xlabel=L"h_I",
        xscale=:log10,
        yscale=:log10,
        ylabel="ϱ",
        legend=:right,
        kwargs...
    )

    h_Is = Base.logrange(1e-4, 1e3, length=100)
    ϱ_theory = zeros(length(h_Is))
    for (k, h_I) in enumerate(h_Is)
        setproperty!(p, :h_I_ini, h_I)
        ϱ_theory[k] = compute_ϱ_ana(p)
    end

    plot!(
        h_Is,
        ϱ_theory;
        label=L"$ϱ_\mathrm{ana}$",
        color=:black,
        linewidth=2,
        kwargs...
    )
    display(current())
    savefig("plots/ice_thickness_convergence.pdf")
end

function plot_a_I_dependence(; iterations=5, kwargs...)
    p = SimulationParameters(Δt_min=600, t_max=3600, Δt_cpl=3600, a_I=1.0, ice_model_type=:temp_feedback)

    a_Is = range(0, 1, 100)
    ϱ_theory = zeros(length(a_Is))
    for (k, a_I) in enumerate(a_Is)
        setproperty!(p, :a_I, a_I)
        ϱ_theory[k] = compute_ϱ_ana(p; s=im * 1e-5)
    end
    plot(
        a_Is,
        ϱ_theory;
        label=L"ϱ_\mathrm{ana}",
        color=:black,
        linewidth=2,
        kwargs...
    )

    a_Is = range(0, 1, 20)
    ϱs_atm = zeros(length(a_Is))
    for (k, a_I) in enumerate(a_Is)
        setproperty!(p, :a_I, a_I)
        _, ϱs_atm[k], _ = run_simulation(p, iterations=iterations)
    end
    plot!(
        a_Is,
        ϱs_atm;
        label=L"$ϱ_\mathrm{num}$",
        markershape=:x,
        color=:black,
        linestyle=:dash,
        linewidth=2,
        legendfontsize=12,
        xlabel=L"a_I",
        ylabel=L"ϱ",
        legend=:bottomright,
        kwargs...
    )

    display(current())
    savefig("plots/ice_a_i_dependence.pdf")
end


function plot_ice_Δt_cpl_convergence(; iterations=10, ice_model_type=:temp_feedback, kwargs...)
    p = SimulationParameters(a_I=1.0, ice_model_type=ice_model_type, Δt_min=10)
    Δt_cpls = Base.logrange(1e1, 1e5, length=5)
    ϱs_atm = zeros(length(Δt_cpls))
    for (k, Δt_cpl) in enumerate(Δt_cpls)
        setproperty!(p, :Δt_cpl, Δt_cpl)
        setproperty!(p, :t_max, Δt_cpl)
        _, ϱs_atm[k], _ = run_simulation(p, iterations=iterations)
    end
    unstable_atm_indices = isinf.(ϱs_atm)
    ϱs_atm[unstable_atm_indices] .= NaN

    plot(
        Δt_cpls,
        ϱs_atm;
        label=L"$ϱ_\mathrm{num}$",
        markershape=:x,
        color=:black,
        linewidth=2,
        legendfontsize=12,
        xlabel=L"Δt_{cpl}",
        xscale=:log10,
        ylabel="ϱ",
        legend=:right,
        kwargs...
    )

    Δt_cpls = Base.logrange(1e1, 1e5, length=100)
    ρ_theory = zeros(length(Δt_cpls))
    for (k, Δt_cpl) in enumerate(Δt_cpls)
        setproperty!(p, :t_max, Δt_cpl)
        ρ_theory[k] = compute_ϱ_ana(p)
    end

    plot!(
        Δt_cpls,
        ρ_theory;
        label=L"$ϱ_\mathrm{ana}$",
        color=:black,
        linewidth=2,
        kwargs...
    )
    display(current())
    savefig("plots/ice_Δt_cpl_convergence.pdf")
end
