using clima_playground
using Plots
using LaTeXStrings
import ClimaCoupler: Interfacer

function plot_C_H_AO_dependence()
    params = SimulationParameters()
    L_AOs = vec(-200:200)
    C_H_AO = similar(L_AOs)
    for (j, L_AO) in enumerate(L_AOs)
        C_H_AO[j] = compute_C_H_AO(params; L_AO=L_AO)
    end
    plot(
        L_AOs,
        C_H_AO,
        xlabel=L"L_{AO}",
        ylabel=L"C_{H,AO}",
        label="",
        color=:black,
    )
end

function plot_C_H_AI_dependence()
    a_Is = 0:0.01:1
    L_AIs = vec(-200:200)
    params = SimulationParameters()
    C_H_AI = zeros(length(a_Is), length(L_AIs))
    for (i, a_I) in enumerate(a_Is)
        params.a_I = a_I
        for (j, L_AI) in enumerate(L_AIs)
            C_H_AI[i, j] = compute_C_H_AI(params; L_AI=L_AI)
        end
    end
    plot()
    surface(L_AIs, a_Is, C_H_AI, xlabel=L"L_{AI}", ylabel=L"a_I", zlabel=L"C_{H,AI}")
end


function plot_ϱ_AO()
    ωs = Base.logrange(1e-6, 1e6, length=100)
    ϱs = similar(ωs)
    params = SimulationParameters()

    for (j, ω) in enumerate(ωs)
        ϱs[j] = compute_ϱ_ana(params; s=im * ω)
    end
    ω_max = ωs[argmax(ϱs)]
    @info "Supremum at ω=$ω_max"

    plot(
        ωs,
        ϱs,
        xlabel="ω",
        ylabel="ϱ(iω)",
        yscale=:log10,
        xscale=:log10,
        color=:black,
        legend=:false,
    )
end


function plot_Δt_cpl_dependence(; plot_title="Δt_cpl_dependence", kwargs...)
    Δt_cpl = [400, 1200, 3600, 10800]
    p = SimulationParameters(Δt_min=400; kwargs...)
    ϱs_atm = similar(Δt_cpl, Float64)
    p.n_t_A = 1
    p.n_t_O = 1

    plot()
    n_z = [(50, 12), (200, 50), (800, 200)]
    for (n_A, n_O) in n_z
        p.n_A = n_A
        p.n_O = n_O
        for (k, var) in enumerate(Δt_cpl)
            p.Δt_min = var
            p.Δt_cpl = var
            p.t_max = var
            _, ϱs_atm[k], _ = run_simulation(p, iterations=10)
        end
        plot!(
            Δt_cpl,
            ϱs_atm,
            label="n_A = $n_A",
            markershape=:x,
            linewidth=2,
        )
    end

    finely_spaced_var = Base.logrange(Δt_cpl[1], Δt_cpl[end], length=100)

    ϱs_analytic = similar(finely_spaced_var)
    for (k, var) in enumerate(finely_spaced_var)
        p.Δt_cpl = var
        p.t_max = var
        ω_min = im * π / p.Δt_cpl
        ω_max = im * π / (p.Δt_min / p.n_t_A)
        ϱs_analytic[k] = max(compute_ϱ_ana(p, s=ω_min), compute_ϱ_ana(p, s=ω_max))
    end
    plot!(
        finely_spaced_var[ϱs_analytic.>0],
        ϱs_analytic[ϱs_analytic.>0],
        label=L"$ϱ_\mathrm{ana}$",
        linewidth=2,
        color=:black,
    )

    plot!(;
        fontsize=18,
        xlabel=L"$\Delta t_{cpl}$",
        ylabel=L"ϱ",
        xscale=:log10,
        yscale=:log10,
        legend=:bottomright,
    )
    display(current())
    savefig("plots/$plot_title.pdf")
end
