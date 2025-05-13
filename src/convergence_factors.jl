import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities

export compute_ϱ_analytical, compute_ϱ_numerical

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


function compute_ϱ_numerical(atmos_vals_list, ocean_vals_list)
    ϱ_A = []
    ϱ_O = []
    pre_bound_error_A = abs.(atmos_vals_list[1] .- atmos_vals_list[end])
    pre_bound_error_O = abs.(ocean_vals_list[1] .- ocean_vals_list[end])
    for i = 2:length(atmos_vals_list)-1
        bound_error_A = abs.(atmos_vals_list[i] .- atmos_vals_list[end])
        bound_error_O = abs.(ocean_vals_list[i] .- ocean_vals_list[end])

        tols_atm = 100 * eps.(max.(abs.(atmos_vals_list[i]), abs.(atmos_vals_list[end])))
        tols_oce = 100 * eps.(max.(abs.(ocean_vals_list[i]), abs.(ocean_vals_list[end])))

        indices_A = findall(
            (pre_bound_error_A[1:end-1] .>= tols_atm[1:end-1]) .&
            (bound_error_A[1:end-1] .>= tols_atm[1:end-1]),
        )
        indices_O = findall(
            (pre_bound_error_O[1:end-1] .>= tols_oce[1:end-1]) .&
            (pre_bound_error_O[1:end-1] .>= tols_oce[1:end-1]),
        )

        ρ_A_value = norm(bound_error_A[indices_A]) / norm(pre_bound_error_A[indices_A])
        ρ_O_value = norm(bound_error_O[indices_O]) / norm(pre_bound_error_O[indices_O])

        push!(ϱ_A, ρ_A_value)
        push!(ϱ_O, ρ_O_value)

        pre_bound_error_A = bound_error_A
        pre_bound_error_O = bound_error_O
    end
    return ϱ_A, ϱ_O
end

