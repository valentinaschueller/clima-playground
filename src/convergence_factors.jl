import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities


function compute_ρ_analytical(p::SimulationParameters; s=nothing)
    if isnothing(s)
        s = im * π / p.t_max
    end
    σ_o = sqrt(s / p.α_o)
    σ_a = sqrt(s / p.α_a)
    ρ = abs((1 - p.a_i)^2 * p.C_AO^2 / (
        (
            p.k_oce *
            σ_o *
            (1 / tanh(σ_o * (p.h_oce - p.z_0numO))) +
            (1 - p.a_i) * p.C_AO +
            p.a_i * p.C_IO
        ) * (
            p.k_atm *
            σ_a *
            (1 / tanh(σ_a * (p.h_atm - p.z_0numA))) +
            (1 - p.a_i) * p.C_AO +
            p.a_i * p.C_AI
        )
    ))
    return ρ
end


function compute_ρ_numerical(atmos_vals_list, ocean_vals_list)
    ρ_atm = []
    ρ_oce = []
    pre_bound_error_atm = abs.(atmos_vals_list[1] .- atmos_vals_list[end])
    pre_bound_error_oce = abs.(ocean_vals_list[1] .- ocean_vals_list[end])
    for i = 2:length(atmos_vals_list)-1
        bound_error_atm = abs.(atmos_vals_list[i] .- atmos_vals_list[end])
        bound_error_oce = abs.(ocean_vals_list[i] .- ocean_vals_list[end])

        tols_atm = 100 * eps.(max.(abs.(atmos_vals_list[i]), abs.(atmos_vals_list[end])))
        tols_oce = 100 * eps.(max.(abs.(ocean_vals_list[i]), abs.(ocean_vals_list[end])))

        indices_atm = findall(
            (pre_bound_error_atm[1:end-1] .>= tols_atm[1:end-1]) .&
            (bound_error_atm[1:end-1] .>= tols_atm[1:end-1]),
        )
        indices_oce = findall(
            (pre_bound_error_oce[1:end-1] .>= tols_oce[1:end-1]) .&
            (pre_bound_error_oce[1:end-1] .>= tols_oce[1:end-1]),
        )

        ρ_atm_value = norm(bound_error_atm[indices_atm]) / norm(pre_bound_error_atm[indices_atm])
        ρ_oce_value = norm(bound_error_oce[indices_oce]) / norm(pre_bound_error_oce[indices_oce])

        push!(ρ_atm, ρ_atm_value)
        push!(ρ_oce, ρ_oce_value)

        pre_bound_error_atm = bound_error_atm
        pre_bound_error_oce = bound_error_oce
    end
    return ρ_atm, ρ_oce
end

"""
Computes the convergence factors based on the parameters in physical_values and for different values of the variable `var`.

**Arguments:**

-`physical_values::Dict`: Can be defined using `define_realistic_vals()`.
-`vars::Array`: Values for the variable.
-`var_name::String`: Name of the variable, for plotting.

**Optional Keyword Arguments:**

-`iterations::Int`: Number of Schwarz iterations, default: 10.
-`a_i_variable::Array`: Values for `a_i` if this should vary as well, default: `nothing`.
-`analytic::Boolean`: Whether to compute the analytical convergence factor, default: `false`.
-`log_scale::Boolean`: Whether the variable is in logarithmic scale, needed to define a range for the corresponding analytical variable, default: `false`.

"""
function get_ρs_one_variable(
    p::SimulationParameters,
    vars,
    var_name;
    iterations=10,
    a_i_variable=nothing,
    analytic=false,
    log_scale=false,
)
    vary_a_i = !isnothing(a_i_variable)
    a_i_variable = vary_a_i ? a_i_variable : [p.a_i]
    ρs_atm = zeros(length(a_i_variable), length(vars))
    ρs_oce = zeros(length(a_i_variable), length(vars))

    if analytic
        if log_scale
            variable2_range = exp10.(range(log10(vars[1]), log10(vars[end]), length=100))
        else
            variable2_range = range(vars[1], vars[end], length=100)
        end
        ρs_analytic = zeros(length(a_i_variable), length(variable2_range))
    else
        variable2_range = nothing
    end

    if var_name == "h_atm"
        Δz_atm =
            (p.h_atm - p.z_0numA) / p.n_atm
    elseif var_name == "h_oce"
        Δz_oce =
            (p.h_oce - p.z_0numO) / p.n_oce
    end

    for (j, a_i) in enumerate(a_i_variable)
        p.a_i = a_i
        for (k, var) in enumerate(vars)
            setproperty!(p, var_name, var)
            restore_physical_values!(p)
            if var_name == :Δt_cpl
                p.t_max = var
            elseif var_name == "h_atm"
                p.n_atm = Int((var - p.z_0numA) / Δz_atm)
            elseif var_name == "h_oce"
                p.n_oce = Int((var - p.z_0numO) / Δz_oce)
            end
            _, ρ_atm, ρ_oce = run_simulation(p, iterations=iterations)
            ρs_atm[j, k], ρs_oce[j, k] = extract_ρ(ρ_atm, ρ_oce)
        end
        if analytic
            for (k, var) in enumerate(variable2_range)
                setproperty!(p, var_name, var)
                restore_physical_values!(p)
                if var_name == :Δt_cpl
                    p.t_max = var
                end
                ρ_analytic = compute_ρ_analytical(p)
                ρs_analytic[j, k] = ρ_analytic
            end
        else
            ρs_analytic = nothing
        end
    end
    return ρs_atm, ρs_oce, variable2_range, ρs_analytic
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
    domain = var1 == :n_atm ? "atm" : "oce"
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
