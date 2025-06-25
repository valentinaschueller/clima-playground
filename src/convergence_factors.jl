import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities

export compute_ϱ_AO, compute_ϱ_numerical, compute_ϱ_AI, compute_ϱ_mixed

function χ_A(p::SimulationParameters, s)
    return tanh((p.h_A - p.z_A0) * sqrt(s / p.α_A))
end

function χ_O(p::SimulationParameters, s)
    return tanh((p.h_O - p.z_O0) * sqrt(s / p.α_O))
end

function compute_ϱ_AO(p::SimulationParameters; s=nothing)
    if isnothing(s)
        s = im * π / p.t_max
    end
    factor = p.k_A / p.k_O * sqrt(p.α_O / p.α_A)
    ν_A = p.k_A / (p.C_AO * sqrt(p.α_A))
    ϱ = factor * abs(χ_A(p, s) / χ_O(p, s)) / abs(1 - ν_A * sqrt(s) * χ_A(p, s))
    return ϱ
end

function compute_ϱ_AI(p::SimulationParameters; s=nothing)
    if isnothing(s)
        s = im * π / p.t_max
    end
    factor = p.C_AI / ((p.k_I / p.h_I_ini) + p.ϵ * p.B + p.C_AI)
    ν_AI = p.k_A / (p.C_AI * sqrt(p.α_A))
    ϱ = factor * abs(1 / (1 - ν_AI * sqrt(s) * χ_A(p, s)))
    return ϱ
end

function compute_ϱ_mixed(p::SimulationParameters; s=nothing)
    if p.ice_model_type == :constant
        ϱ_AI = 0.0
    else
        ϱ_AI = compute_ϱ_AI(p; s)
    end
    ϱ_AO = compute_ϱ_AO(p; s)
    ϱ_mixed = p.a_I * ϱ_AI + (1 - p.a_I) * ϱ_AO
    return ϱ_mixed
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

