using clima_playground
using Plots
using LaTeXStrings


function plot_aspect_ratio_dependence(; plot_title="aspect_ratio", kwargs...)
    p = SimulationParameters{Float64}(Δt_min=600, n_t_A=10, t_max=3600, Δt_cpl=3600, ice_model_type=:constant; kwargs...)
    
    plot()
    n_z = [(200, 400), (200, 200), (200, 100), (200, 50), (400, 50), (800, 50), (1600, 50), (1600, 25)]
    r = [0.25 * n_A / n_O for (n_A, n_O) in n_z]
    ϱs_atm = similar(r, Float64)
    k = 1
    for (n_A, n_O) in n_z
        p.n_A = n_A
        p.n_O = n_O
        _, ϱs_atm[k], _ = run_simulation(p, iterations=10)
        k += 1
    end
    plot!(
            r,
            ϱs_atm,
            label=L"$ϱ_\mathrm{num}$",
            markershape=:x,
            linewidth=2,
            color=:black,
        )

    ω_min = im * π / p.Δt_cpl
    ω_max = im * π / (p.Δt_min / p.n_t_A)
    ϱ_analytic = max(compute_ϱ_ana(p, s=ω_min), compute_ϱ_ana(p, s=ω_max))
    hline!(
        [ϱ_analytic],
        label=L"$ϱ_\mathrm{ana}$",
        linewidth=2,
        color=:black,
    )

    plot!(;
        fontsize=18,
        xlabel="Aspect Ratio",
        ylabel=L"ϱ",
        yscale=:log10,
        xscale=:log2,
        legend=:bottomright,
        color=:black,
        ylim=[1e-5, 5e-3],
    )
    display(current())
    savefig("plots/$plot_title.pdf")
end

function plot_aspect_ratio_dependence_same_dt(; plot_title="aspect_ratio_same_dt", kwargs...)
    p = SimulationParameters{Float64}(Δt_min=600, n_t_A=1, n_t_O=1, t_max=600, Δt_cpl=600, ice_model_type=:constant; kwargs...)
    
    plot()
    n_z = [(200, 400), (200, 200), (200, 100), (200, 50), (400, 50), (800, 50), (1600, 50), (1600, 25)]
    r = [0.25 * n_A / n_O for (n_A, n_O) in n_z]
    ϱs_atm = similar(r, Float64)
    k = 1
    for (n_A, n_O) in n_z
        p.n_A = n_A
        p.n_O = n_O
        _, ϱs_atm[k], _ = run_simulation(p, iterations=10)
        k += 1
    end
    plot!(
            r,
            ϱs_atm,
            label=L"$ϱ_\mathrm{num}$",
            markershape=:x,
            linewidth=2,
            color=:black,
        )

    ω_min = im * π / p.Δt_cpl
    ω_max = im * π / (p.Δt_min / p.n_t_A)
    ϱ_analytic = max(compute_ϱ_ana(p, s=ω_min), compute_ϱ_ana(p, s=ω_max))
    hline!(
        [ϱ_analytic],
        label=L"$ϱ_\mathrm{ana}$",
        linewidth=2,
        color=:black,
    )

    plot!(;
        fontsize=18,
        xlabel="Aspect Ratio",
        ylabel=L"ϱ",
        yscale=:log10,
        xscale=:log2,
        legend=:bottomright,
        ylim=[1e-5, 5e-3],
    )
    display(current())
    savefig("plots/$plot_title.pdf")
end

