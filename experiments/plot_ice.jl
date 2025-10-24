using clima_playground
using Plots
using LaTeXStrings
import ClimaCore as CC

function plot_ice_seb_results()
    params = SimulationParameters{Float64}(C_H_AI=1.4e-3, h_I_ini=1.0)
    T_As = vec(range(270, 290, length=100))
    T_ice = similar(T_As)
    for (i, T_A) in enumerate(T_As)
        params.T_A = T_A
        T_ice[i] = T_Is(params)
    end
    plot(T_As, T_ice, color=:black, xlabel=L"T_A", ylabel=L"T_{I,s}", legend=false)
    display(current())
end

function plot_ice_thickness_convergence(; plot_title="ice_thickness_convergence", iterations=5, kwargs...)
    p = SimulationParameters{Float64}(Δt_min=600, t_max=3600, Δt_cpl=3600, n_t_A=10, a_I=1.0, C_AI=1.82, ice_model_type=:thickness_feedback)
    h_Is = Base.logrange(5e-4, 5e1, length=10)
    ϱs_atm = similar(h_Is)
    for (k, h_I) in enumerate(h_Is)
        setproperty!(p, :h_I_ini, h_I)
        _, ϱs_atm[k], _ = run_simulation(p, iterations=iterations)
    end

    p.ice_model_type = :temp_feedback
    ϱs_atm_lin = similar(h_Is)
    for (k, h_I) in enumerate(h_Is)
        setproperty!(p, :h_I_ini, h_I)
        _, ϱs_atm_lin[k], _ = run_simulation(p, iterations=iterations)
    end

    finely_spaced_var = Base.logrange(h_Is[1], h_Is[end], length=100)
    ϱ_theory = similar(finely_spaced_var)
    for (k, h_I) in enumerate(finely_spaced_var)
        setproperty!(p, :h_I_ini, h_I)
        ϱ_theory[k] = compute_ϱ_ana(p)
    end

    pgfplotsx()
    plot(
        finely_spaced_var,
        ϱ_theory;
        label=L"$\varrho_{AI}(\omega_\mathrm{max})$",
        color=:black,
    )

    plot!(
        h_Is,
        ϱs_atm;
        label=L"\varrho_\mathrm{num}, h_I=h_I(t)",
        markershape=:circle,
        ms=4,
        ls=:dot,
        color=:black,
    )

    plot!(
        h_Is,
        ϱs_atm_lin;
        label=L"\varrho_\mathrm{num}, h_I=h_0",
        markershape=:x,
        color=:black,
        linestyle=:dash,
        ms=6,
    )
    plot!(
        labelfontsize=20,
        tickfontsize=20,
        legendfontsize=18,
        xlabel=L"h_0",
        xscale=:log10,
        yscale=:log10,
        ylabel=L"\varrho",
        legend=:right,
        ylim=[6e-5, 1.5],
        yticks=[1e-4, 1e-2, 1],
        yminorticks=1,
        xticks=[1e-3, 1e-2, 1e-1, 1, 1e1],
        kwargs...
    )

    display(current())
    savefig("plots/$plot_title.tikz")
end

function plot_a_I_dependence(; plot_title="ice_a_i_dependence", iterations=5, kwargs...)
    p = SimulationParameters{Float64}(Δt_min=600, t_max=3600, Δt_cpl=3600, a_I=1.0, n_t_A=10, C_AI=1.82, C_AO=1.3, ice_model_type=:thickness_feedback)

    a_Is = range(0, 1, 500)
    ϱ_theory = similar(a_Is)
    for (k, a_I) in enumerate(a_Is)
        setproperty!(p, :a_I, a_I)
        ϱ_theory[k] = compute_ϱ_ana(p)
    end

    pgfplotsx()
    plot(
        a_Is,
        ϱ_theory;
        label=L"\varrho(\omega_\mathrm{max})",
        color=:black,
        kwargs...
    )

    a_Is = range(0, 1, 15)
    ϱs_atm = similar(a_Is)
    for (k, a_I) in enumerate(a_Is)
        setproperty!(p, :a_I, a_I)
        _, ϱs_atm[k], _ = run_simulation(p, iterations=iterations)
    end
    plot!(
        a_Is,
        ϱs_atm;
        label=L"$\varrho_\mathrm{num}$",
        markershape=:x,
        ms=6,
        color=:black,
        labelfontsize=22,
        tickfontsize=22,
        legendfontsize=18,
        xlabel=L"a_I",
        ylabel=L"\varrho",
        legend=:right,
        yscale=:log10,
        ylim=[1e-6, 1],
        yminorticks=1,
        yticks=[1e-6, 1e-4, 1e-2, 1],
        kwargs...
    )

    display(current())
    savefig("plots/$plot_title.tikz")
end

function plot_C_AX_dependence(; plot_title="C_AX_dependence", kwargs...)
    bulk_coeffs = Base.logrange(1e-2, 1e2, length=10)
    finely_spaced_var = Base.logrange(bulk_coeffs[1], bulk_coeffs[end], length=100)
    p = SimulationParameters{Float64}(Δt_min=600, n_t_A=10, t_max=3600, Δt_cpl=3600, a_I=1.0, ice_model_type=:thickness_feedback; kwargs...)
    ϱs_atm = similar(bulk_coeffs)

    for (k, C_AI) in enumerate(bulk_coeffs)
        p.C_AI = C_AI
        _, ϱs_atm[k], _ = run_simulation(p, iterations=5)
    end

    ϱs_analytic = similar(finely_spaced_var)
    for (k, C_AI) in enumerate(finely_spaced_var)
        p.C_AI = C_AI
        ϱs_analytic[k] = compute_ϱ_ana(p)
    end

    pgfplotsx()
    plot(
        finely_spaced_var,
        ϱs_analytic,
        label=L"$\varrho_\mathrm{AI}(\omega_\mathrm{max})$",
        color=:black,
        linestyle=:dash,
    )
    plot!(
        bulk_coeffs,
        ϱs_atm,
        label=L"\varrho_\mathrm{num}, a_I=1",
        color=:black,
        markershape=:x,
        ms=6,
        linestyle=:dash,
    )

    p = SimulationParameters{Float64}(Δt_min=10, t_max=1000, Δt_cpl=1000, a_I=0.0, ice_model_type=:temp_feedback; kwargs...)
    ϱs_atm = similar(bulk_coeffs)

    for (k, C_AO) in enumerate(bulk_coeffs)
        p.C_AO = C_AO
        _, ϱs_atm[k], _ = run_simulation(p, iterations=5)
    end

    ϱs_analytic = similar(finely_spaced_var)
    for (k, C_AO) in enumerate(finely_spaced_var)
        p.C_AO = C_AO
        ϱs_analytic[k] = compute_ϱ_ana(p)
    end
    plot!(
        finely_spaced_var,
        ϱs_analytic,
        label=L"$\varrho_\mathrm{AO}(\omega_\mathrm{max})$",
        color=:black,
    )
    plot!(
        bulk_coeffs,
        ϱs_atm,
        label=L"\varrho_\mathrm{num}, a_I=0",
        color=:black,
        markershape=:x,
        ms=6,
    )
    plot!(;
        labelfontsize=22,
        tickfontsize=22,
        legendfontsize=18,
        legendcolumns=2,
        xlabel=L"$C_{AO}, C_{AI}$",
        ylabel=L"\varrho",
        yscale=:log10,
        xscale=:log10,
        legend=:right,
        yticks=[1e-6, 1e-4, 1e-2, 1],
        ylim=[5e-7, 1.5],
    )
    display(current())
    savefig("plots/$plot_title.tikz")
end
