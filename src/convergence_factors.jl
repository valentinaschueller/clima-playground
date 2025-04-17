import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities

"""
Computes the analytical convergence factor based on the parameters in physical_values.

**Arguments:**

-`physical_values::Dict`: Can be defined using `define_realistic_vals()`.

"""
function compute_ρ_analytical(physical_values)
    σ_o = sqrt(abs(physical_values[:w_min]) / (2 * physical_values[:α_o])) * (1 + im)
    σ_a = sqrt(abs(physical_values[:w_min]) / (2 * physical_values[:α_a])) * (1 + im)
    η_AO =
        physical_values[:C_AO] *
        abs(physical_values[:u_atm] - physical_values[:u_oce]) *
        physical_values[:ρ_atm] *
        physical_values[:c_atm]
    η_OI =
        physical_values[:C_OI] *
        abs(physical_values[:u_oce]) *
        physical_values[:ρ_oce] *
        physical_values[:c_oce]
    η_AI =
        physical_values[:C_AI] *
        abs(physical_values[:u_atm]) *
        physical_values[:ρ_atm] *
        physical_values[:c_atm]
    return abs(
        (1 - physical_values[:a_i])^2 * η_AO^2 / (
            (
                physical_values[:k_oce] *
                σ_o *
                (1 / tanh(σ_o * (physical_values[:h_oce] - physical_values[:z_0numO]))) +
                (1 - physical_values[:a_i]) * η_AO +
                physical_values[:a_i] * η_OI
            ) * (
                physical_values[:k_atm] *
                σ_a *
                (1 / tanh(σ_a * (physical_values[:h_atm] - physical_values[:z_0numA]))) +
                (1 - physical_values[:a_i]) * η_AO +
                physical_values[:a_i] * η_AI
            )
        ),
    )
end

"""
Computes the numerical convergence factor based on the parameters in physical_values.

**Arguments:**

-`physical_values::Dict`: Can be defined using `define_realistic_vals()`.

**Optional Keyword Arguments:**

-`iterations::Int`: Number of Schwarz iterations, default: 1.

"""
function compute_ρ_numerical(
    physical_values;
    iterations=1,
)
    cs = get_coupled_sim(physical_values)
    ρ_atm, ρ_oce = solve_coupler!(
        cs,
        iterations=iterations,
    )
    return ρ_atm, ρ_oce
end

"""
Computes the convergence factors based on the parameters in physical_values and for different values of the variable `var`.

**Arguments:**

-`physical_values::Dict`: Can be defined using `define_realistic_vals()`.
-`vars::Array`: Values for the variable.
-`var_name::String`: Name of the variable, for plotting.

**Optional Keyword Arguments:**

-`iterations::Int`: Number of Schwarz iterations, default: 1.
-`a_i_variable::Array`: Values for `a_i` if this should vary as well, default: `nothing`.
-`analytic::Boolean`: Whether to compute the analytical convergence factor, default: `false`.
-`log_scale::Boolean`: Whether the variable is in logarithmic scale, needed to define a range for the corresponding analytical variable, default: `false`.

"""
function get_ρs_one_variable(
    physical_values,
    vars,
    var_name;
    iterations=1,
    a_i_variable=nothing,
    analytic=false,
    log_scale=false,
)
    vary_a_i = !isnothing(a_i_variable)
    a_i_variable = vary_a_i ? a_i_variable : [physical_values[:a_i]]
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
            (physical_values[:h_atm] - physical_values[:z_0numA]) / physical_values[:n_atm]
    elseif var_name == "h_oce"
        Δz_oce =
            (physical_values[:h_oce] - physical_values[:z_0numO]) / physical_values[:n_oce]
    end

    for (j, a_i) in enumerate(a_i_variable)
        physical_values[:a_i] = a_i
        correct_for_a_i!(physical_values)
        for (k, var) in enumerate(vars)
            physical_values[Symbol(var_name)] = var
            if var_name == "a_i"
                correct_for_a_i!(physical_values)
            elseif var_name == "h_atm"
                physical_values[:n_atm] = Int((var - physical_values[:z_0numA]) / Δz_atm)
            elseif var_name == "h_oce"
                physical_values[:n_oce] = Int((var - physical_values[:z_0numO]) / Δz_oce)
            elseif var_name == "t_max"
                physical_values[:Δt_cpl] = var
            end
            ρ_atm, ρ_oce =
                compute_ρ_numerical(physical_values, iterations=iterations)
            ρs_atm[j, k], ρs_oce[j, k] =
                extract_ρ(ρ_atm, ρ_oce)
        end
        if analytic
            for (k, var) in enumerate(variable2_range)
                physical_values[Symbol(var_name)] = var
                if var_name == "a_i"
                    correct_for_a_i!(physical_values)
                elseif var_name == "t_max"
                    physical_values[:w_min] = π / var
                end
                ρ_analytic = compute_ρ_analytical(physical_values)
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
function stability_check(physical_values, n_zs, n_ts, var1_name, var2_name)
    domain = var1_name[end-2:end]
    unstable_matrix = zeros(length(n_ts), length(n_zs))

    for (i, n_z) in enumerate(n_zs)
        for (j, n_t) in enumerate(n_ts)
            physical_values[Symbol(var1_name)] = n_z
            physical_values[Symbol(var2_name)] = n_t

            cs = get_coupled_sim(physical_values)

            lower_limit_temp, upper_limit_temp = initial_value_range(cs)

            if domain == "atm"
                Interfacer.step!(cs.model_sims.atmos_sim, physical_values[:Δt_cpl])
                states = copy(cs.model_sims.atmos_sim.integrator.sol.u)
            elseif domain == "oce"
                Interfacer.step!(cs.model_sims.ocean_sim, physical_values[:Δt_cpl])
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
