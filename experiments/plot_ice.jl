using clima_playground
using Plots
using LaTeXStrings
import ClimaCore as CC

function plot_ice_seb_results()
    params = SimulationParameters(C_H_AI=1.4e-3)
    T_As = vec(range(270, 290, length=100))
    T_ice = similar(T_As)
    for (i, T_A) in enumerate(T_As)
        params.T_A = T_A
        T_ice[i] = compute_T_Is(params, 1.0)
    end
    plot(T_As, T_ice, color=:black, xlabel=L"T_A", ylabel=L"T_{I,s}", legend=false)
    display(current())
end

function plot_ice_thickness_convergence(; iterations=5, kwargs...)
    p = SimulationParameters(Δt_min=200, t_max=3600, Δt_cpl=3600, a_I=1.0, ice_model_type=:thickness_feedback)
    h_Is = Base.logrange(1e-4, 1e3, length=10)
    ϱs_atm = similar(h_Is)
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
    ϱs_atm = similar(h_Is)
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
    ϱ_theory = similar(h_Is)
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
    ϱ_theory = similar(a_Is)
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
    ϱs_atm = similar(a_Is)
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

function plot_a_I_zoom(; iterations=5, kwargs...)
    p = SimulationParameters(Δt_min=600, t_max=1200, Δt_cpl=1200, a_I=1.0, ice_model_type=:temp_feedback)

    a_Is = range(1e-4, 3e-4, 500)
    ϱs_atm = similar(a_Is)
    for (k, a_I) in enumerate(a_Is)
        setproperty!(p, :a_I, a_I)
        _, ϱs_atm[k], _ = run_simulation(p, iterations=iterations)
    end
    plot(
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
    savefig("plots/ice_a_i_zoom.pdf")
end


function plot_ice_Δt_cpl_convergence(; iterations=10, ice_model_type=:temp_feedback, kwargs...)
    p = SimulationParameters(a_I=1.0, ice_model_type=ice_model_type, Δt_min=10)
    Δt_cpls = Base.logrange(1e1, 1e5, length=5)
    ϱs_atm = similar(Δt_cpls)
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
    ρ_theory = similar(Δt_cpls)
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

function plot_C_AI_dependence(plot_title="C_AI_dependence"; kwargs...)
    C_AIs = Base.logrange(1e-2, 1e2, length=15)
    p = SimulationParameters(Δt_min=10, t_max=1000, Δt_cpl=1000, a_I=1.0, ice_model_type=:temp_feedback; kwargs...)
    ϱs_atm = similar(C_AIs)

    for (k, var) in enumerate(C_AIs)
        p.C_AI = var
        _, ϱs_atm[k], _ = run_simulation(p, iterations=6)
    end

    finely_spaced_var = Base.logrange(C_AIs[1], C_AIs[end], length=100)
    ϱs_analytic = similar(finely_spaced_var)
    for (k, var) in enumerate(finely_spaced_var)
        p.C_AI = var
        ω_min = im * π / p.Δt_cpl
        ω_max = im * π / (p.Δt_min / p.n_t_A)
        ϱs_analytic[k] = max(compute_ϱ_ana(p, s=ω_min), compute_ϱ_ana(p, s=ω_max))
    end
    plot(
        finely_spaced_var,
        ϱs_analytic,
        label=L"$ϱ_\mathrm{ana}$",
        linewidth=2,
        color=:black,
    )
    plot!(
        C_AIs,
        ϱs_atm,
        label=L"$ϱ_\mathrm{num}$",
        color=:black,
        markershape=:x,
        linewidth=2,
    )
    plot!(;
        legendfontsize=12,
        xlabel=L"C_{AI}",
        ylabel="ϱ",
        xscale=:log10,
        yscale=:log10,
        legend=:bottomright,
    )
    display(current())
    savefig("plots/$plot_title.pdf")
end
