"""
Extracts an Array from an Array of FieldVectors.

**Arguments:**

-`field_vecs: Array`: An array of fieldvectors containing temperature values.
-`domain: String`: Either `"atm"` or `"oce"`.

"""
function extract_matrix(field_vecs, domain)
    matrix = []
    for field_vec in field_vecs
        field = (domain == "atm") ? field_vec.atm : field_vec.oce
        values = parent(field)
        push!(matrix, values)
    end
    return hcat(matrix...)
end

function initial_value_range(cs)
    max_value = floatmin()
    min_value = floatmax()
    for model_sim in cs.model_sims
        initial_values = parent(model_sim.Y_init)
        min_value = min(min_value, minimum(initial_values))
        max_value = max(max_value, maximum(initial_values))
    end
    return min_value, max_value
end

"""
Checks if the computed temperatures are reasonable or if the model has gone unstable.

**Arguments:**

-`values: Array`: Array of temperature values.
-`upper_limit: Float64`: Maximum reasonable temperature.
-`lower_limit: Float64`: Minimum reasonable temperature.
"""
function is_stable(values, upper_limit, lower_limit)
    if (
        any(isnan, values) ||
        maximum(values) > upper_limit ||
        minimum(values) < lower_limit
    )
        return false
    end
    return true
end

"""
Checks if the Schwarz iteration should terminated.

**Arguments:**
    
-`bound_atmos_vals: Array`: Atmosphere boundary temperatures.
-`pre_bound_atmos_vals: Array`: Previous iteration atmosphere boundary temperatures.
-`bound_ocean_vals: Array`: Ocean boundary temperatures.
-`pre_bound_ocean_vals: Array`: Previous iteration ocean boundary temperatures.    
-`iter: Int`: Current iteration.

"""
function has_converged(
    bound_atmos_vals,
    pre_bound_atmos_vals,
    bound_ocean_vals,
    pre_bound_ocean_vals
)
    if isnothing(pre_bound_atmos_vals)
        return false
    end
    bound_errors_atm_iter = abs.(bound_atmos_vals .- pre_bound_atmos_vals)
    bound_errors_oce_iter = abs.(bound_ocean_vals .- pre_bound_ocean_vals)
    tols_atm = 100 * eps.(max.(abs.(bound_atmos_vals), abs.(pre_bound_atmos_vals)))
    tols_oce = 100 * eps.(max.(abs.(bound_ocean_vals), abs.(pre_bound_ocean_vals)))
    if all(bound_errors_atm_iter .< tols_atm)
        return true
    elseif all(bound_errors_oce_iter .< tols_oce)
        return true
    else
        return false
    end
end

"""
Special treatment when varying some variables.

**Arguments:**

-`var::Array`: Values for the variable.
-`var_name::String`: The variable name.
-`ρs_oce::Array`: Ocean convergence factors.
-`ρs_atm::Array`: Atmosphere convergence factors.
-`physical_values::Dict`: Can be defined using `define_realistic_vals()`.

**Optional Keyword Arguments:**

-`dims::Int`: Dimensions over which to  reverse the convergence factor matrices, default: `1`.
-`param_analytic::Array`: The analytical values for the variable, default: `nothing`.
-`ρs_analytic::Array`: The analytical convergence factors, default: `nothing`.

"""
function handle_variable(
    var,
    var_name,
    ρs_oce,
    ρs_atm,
    physical_values;
    dims=1,
    param_analytic=nothing,
    ρs_analytic=nothing,
)
    # Special treatment of n_atm and n_oce.
    if var_name == "n_atm"
        var = (physical_values[:h_atm] - physical_values[:z_0numA]) ./ reverse(var)
        if !isnothing(ρs_analytic) && !isnothing(param_analytic)
            param_analytic = reverse(
                (physical_values[:h_atm] - physical_values[:z_0numA]) ./ param_analytic,
            )
        end
    elseif var_name == "n_oce"
        var = (physical_values[:h_oce] - physical_values[:z_0numO]) ./ reverse(var)
        if !isnothing(ρs_analytic) && !isnothing(param_analytic)
            param_analytic = reverse(
                (physical_values[:h_oce] - physical_values[:z_0numO]) ./ param_analytic,
            )
        end
    elseif var_name == "n_t_atm" || var_name == "n_t_oce"
        var = physical_values[:Δt_min] ./ reverse(var)
        if !isnothing(ρs_analytic) && !isnothing(param_analytic)
            param_analytic = reverse(physical_values[:Δt_min] ./ param_analytic)
        end
    end
    if var_name in ["n_atm", "n_oce", "n_t_atm", "n_t_oce"]
        ρs_oce =
            !isnothing(ρs_oce) ? reverse(ρs_oce, dims=dims) : nothing
        ρs_atm =
            !isnothing(ρs_atm) ? reverse(ρs_atm, dims=dims) : nothing
        ρs_analytic =
            !isnothing(ρs_analytic) ? reverse(ρs_analytic, dims=dims) :
            nothing
    end
    return var, ρs_oce, ρs_atm, param_analytic, ρs_analytic
end

""" 
Extracts the first convergence factors.

**Arguments:**

-`ρ_atm::Array`: Atmosphere convergence factor.
-`ρ_oce::Array`: Ocean convergence factor.

"""
function extract_ρ(ρ_atm, ρ_oce)
    if ρ_atm isa AbstractArray
        ρ_atm = !isempty(ρ_atm) ? ρ_atm[1] : NaN
    end
    if ρ_oce isa AbstractArray
        ρ_oce = !isempty(ρ_oce) ? ρ_oce[1] : NaN
    end
    return ρ_atm, ρ_oce
end

function update_ρ_parallel(ρ_atm, ρ_oce)
    ρ_atm = ρ_atm[1:end-1] .* ρ_atm[2:end]
    ρ_oce = ρ_oce[1:end-1] .* ρ_oce[2:end]
    ρ_atm[1:2:end] .= NaN
    ρ_oce[2:2:end] .= NaN
    return ρ_atm, ρ_oce
end
