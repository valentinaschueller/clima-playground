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

function plot_ice_thickness_convergence(; iterations=10, kwargs...)
    p = SimulationParameters(t_max=1000, Δt_cpl=1000, a_I=1.0, ice_model_type=:temp_feedback)
    h_Is = Base.logrange(1e-4, 1e3, length=7)
    ϱs_atm = zeros(length(h_Is))
    for (k, h_I) in enumerate(h_Is)
        setproperty!(p, :h_I_ini, h_I)
        _, ϱ_atm, ϱ_oce = run_simulation(p, iterations=iterations)
        ϱs_atm[k], _ = extract_ρ(ϱ_atm, ϱ_oce)
    end
    unstable_atm_indices = isinf.(ϱs_atm)
    ϱs_atm[unstable_atm_indices] .= NaN

    plot(
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

    h_Is = Base.logrange(1e-4, 1e3, length=100)
    ρ_theory = zeros(length(h_Is))
    for (k, h_I) in enumerate(h_Is)
        setproperty!(p, :h_I_ini, h_I)
        ρ_theory[k] = compute_ϱ_ice(p)
    end

    plot!(
        h_Is,
        ρ_theory;
        label=L"$ϱ_\mathrm{ana}$",
        color=:black,
        linewidth=2,
        kwargs...
    )
    savefig("plots/ice_thickness_convergence.pdf")
end

function plot_a_I_dependence(; iterations=10, kwargs...)
    p = SimulationParameters(t_max=1000, Δt_cpl=1000, a_I=1.0, ice_model_type=:temp_feedback, C_H_IO=7.5e-3)
    a_Is = range(0, 1, 30)
    ϱs_atm = zeros(length(a_Is))
    for (k, a_I) in enumerate(a_Is)
        setproperty!(p, :a_I, a_I)
        _, ϱ_atm, ϱ_oce = run_simulation(p, iterations=iterations)
        ϱs_atm[k], _ = extract_ρ(ϱ_atm, ϱ_oce)
    end

    unstable_atm_indices = isinf.(ϱs_atm)
    ϱs_atm[unstable_atm_indices] .= NaN

    plot(
        a_Is,
        ϱs_atm;
        label=L"$ϱ_\mathrm{num}$",
        markershape=:x,
        color=:black,
        linewidth=2,
        legendfontsize=12,
        xlabel=L"a_I",
        yscale=:log10,
        ylabel="ϱ",
        legend=:right,
        kwargs...
    )
    savefig("plots/ice_a_i_dependence.pdf")
end


function plot_ice_Δt_cpl_convergence(; iterations=10, kwargs...)
    p = SimulationParameters(a_I=1.0, ice_model_type=:temp_feedback)
    Δt_cpls = Base.logrange(1e1, 1e5, length=5)
    ϱs_atm = zeros(length(Δt_cpls))
    for (k, Δt_cpl) in enumerate(Δt_cpls)
        setproperty!(p, :Δt_cpl, Δt_cpl)
        setproperty!(p, :t_max, Δt_cpl)
        _, ϱ_atm, ϱ_oce = run_simulation(p, iterations=iterations)
        ϱs_atm[k], _ = extract_ρ(ϱ_atm, ϱ_oce)
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
        ρ_theory[k] = compute_ϱ_ice(p)
    end

    plot!(
        Δt_cpls,
        ρ_theory;
        label=L"$ϱ_\mathrm{ana}$",
        color=:black,
        linewidth=2,
        kwargs...
    )
    savefig("plots/ice_Δt_cpl_convergence.pdf")
end
