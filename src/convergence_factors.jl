import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities

export compute_ϱ_numerical, compute_ϱ_ana

function χ_A(p::SimulationParameters, s)
    return tanh((p.h_A - p.z_A0) * sqrt(s / p.α_A))
end

function χ_O(p::SimulationParameters, s)
    return tanh((p.h_O - p.z_O0) * sqrt(s / p.α_O))
end

function compute_ϱ_ana(p::SimulationParameters; s=nothing)
    if isnothing(s)
        s = im * π / p.t_max
    end
    ζ = p.C_AI / (p.k_I / p.h_I_ini + p.ϵ * p.B + p.C_AI)
    frac_in_num = (p.k_A / sqrt(p.α_A)) * sqrt(s) * χ_A(p, s) / (p.a_I * p.C_IO - (p.k_O / sqrt(p.α_O)) * sqrt(s) * χ_O(p, s))
    num = p.a_I * p.C_AI * ζ + (1 - p.a_I)^2 * p.C_AO * frac_in_num
    den = (p.k_A / sqrt(p.α_A)) * sqrt(s) * χ_A(p, s) - p.a_I * p.C_AI - (1 - p.a_I) * p.C_AO
    return abs(num / den)
end

function compute_ϱ_numerical(coupling_variable)
    ϱ = []
    e_old = abs.(coupling_variable[1] .- coupling_variable[end])
    for i = 2:length(coupling_variable)-1
        e_new = abs.(coupling_variable[i] .- coupling_variable[end])
        ϱ_i = norm(e_new) / norm(e_old)
        push!(ϱ, ϱ_i)

        e_old = e_new
    end
    return mean_ϱ(ϱ)
end

