import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities


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



"""
Computes for which `Δzᴬ` and `Δtᴬ` or `Δzᴼ` and `Δtᴼ` the model becomes unstable.

**Arguments:**

-`physical_values: Dict`: Can be defined using `define_realistic_vals()`.
-`n_zs::Array`: Values for the number of spatial gridpoints.
-`n_ts::Array`: Values for the number of timepoints.
-`var1_name::String`: Either `"Δzᴬ"` or `"Δzᴼ"`.
-`var1_name::String`: Either `"Δtᴬ"` or `"Δtᴼ"`.

"""
function stability_check(p::SimulationParameters, n_zs, n_ts, var1, var2)
    domain = var1 == :n_A ? "atm" : "oce"
    unstable_matrix = zeros(length(n_ts), length(n_zs))

    for (i, n_z) in enumerate(n_zs)
        for (j, n_t) in enumerate(n_ts)
            setproperty!(p, var1, n_z)
            setproperty!(p, var2, n_t)

            cs = get_coupled_sim(p)

            lower_limit_temp, upper_limit_temp = initial_value_range(cs)

            if domain == "atm"
                Interfacer.step!(cs.model_sims.atmos_sim, p.Δt_cpl)
                states = copy(cs.model_sims.atmos_sim.integrator.sol.u)
            else
                Interfacer.step!(cs.model_sims.ocean_sim, p.Δt_cpl)
                states = copy(cs.model_sims.ocean_sim.integrator.sol.u)
            end
            vals = extract_matrix(states, domain)
            if !is_stable(vals, upper_limit_temp, lower_limit_temp)
                unstable_matrix[i, j] = Inf
            else
                unstable_matrix[i, j] = NaN
            end
        end
    end
    return unstable_matrix
end
