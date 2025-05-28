import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities

export compute_ϱ_analytical, compute_ϱ_numerical, compute_ϱ_ice

function compute_ϱ_analytical(p::SimulationParameters; s=nothing)
    if isnothing(s)
        s = im * π / p.t_max
    end
    σ_o = sqrt(s / p.α_O)
    σ_a = sqrt(s / p.α_A)
    ϱ = abs((1 - p.a_I)^2 * p.C_AO^2 / (
        (
            p.k_O *
            σ_o *
            (1 / tanh(σ_o * (p.h_O - p.z_O0))) +
            (1 - p.a_I) * p.C_AO +
            p.a_I * p.C_IO
        ) * (
            p.k_A *
            σ_a *
            (1 / tanh(σ_a * (p.h_A - p.z_A0))) +
            (1 - p.a_I) * p.C_AO +
            p.a_I * p.C_AI
        )
    ))
    return ϱ
end

function compute_ϱ_ice(p::SimulationParameters; s=nothing)
    if isnothing(s)
        s = im * π / p.t_max
    end
    factor = p.C_AI / ((p.k_I / p.h_I_ini) + p.ϵ * p.B + p.C_AI)
    ξ_A = (p.h_A - p.z_A0) / sqrt(p.α_A)
    ν_AI = p.k_A / (p.C_AI * sqrt(p.α_A))
    ϱ = factor * abs(1 / (1 - ν_AI * sqrt(s) * tanh(ξ_A * sqrt(s))))
    return ϱ
end

function compute_ϱ_numerical(coupling_variable)
    ϱ = []
    e_old = abs.(coupling_variable[1] .- coupling_variable[end])
    for i = 2:length(coupling_variable)-1
        e_new = abs.(coupling_variable[i] .- coupling_variable[end])

        tols = 100 * eps.(max.(abs.(coupling_variable[i]), abs.(coupling_variable[end])))

        indices = findall(
            (e_old[1:end-1] .>= tols[1:end-1]) .&
            (e_new[1:end-1] .>= tols[1:end-1]),
        )

        ϱ_i = norm(e_new[indices]) / norm(e_old[indices])
        push!(ϱ, ϱ_i)

        e_old = e_new
    end
    return mean_ϱ(ϱ)
end

