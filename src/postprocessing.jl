"""
Extracts an Array from an Array of FieldVectors.

**Arguments:**

-`field_vecs: Array`: An array of fieldvectors containing temperature values.
-`domain: String`: Either `"atm"` or `"oce"`.

"""
function extract_matrix(field_vecs, domain)
    matrix = []
    for field_vec in field_vecs
        field = field_vec.data
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
Extracts the first convergence factors.

**Arguments:**

-`ρ_A::Array`: Atmosphere convergence factor.
-`ρ_O::Array`: Ocean convergence factor.

"""
function extract_ρ(ρ_A, ρ_O)
    if ρ_A isa AbstractArray
        ρ_A = !isempty(ρ_A) ? ρ_A[1] : NaN
    end
    if ρ_O isa AbstractArray
        ρ_O = !isempty(ρ_O) ? ρ_O[1] : NaN
    end
    return ρ_A, ρ_O
end

function update_ρ_parallel(ρ_A, ρ_O)
    ρ_A = ρ_A[1:end-1] .* ρ_A[2:end]
    ρ_O = ρ_O[1:end-1] .* ρ_O[2:end]
    ρ_A[1:2:end] .= NaN
    ρ_O[2:2:end] .= NaN
    return ρ_A, ρ_O
end
