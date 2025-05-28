export check_stability, has_converged, mean_ϱ, update_ρ_parallel, UnstableError

struct UnstableError <: Exception
end

function check_stability(values, value_range=nothing)
    if any(isnan, values)
        throw(UnstableError())
    end
    if !isnothing(value_range)
        if (
            maximum(values) > maximum(value_range) ||
            minimum(values) < minimum(value_range)
        )
            throw(UnstableError())
        end
    end
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
    if all(bound_errors_atm_iter .< tols_atm) && all(bound_errors_oce_iter .< tols_oce)
        return true
    else
        return false
    end
end

function mean_ϱ(ϱ_array)
    non_nan_values = filter(!isnan, ϱ_array)
    if isempty(non_nan_values)
        return NaN
    end
    return mean(non_nan_values)
end

function update_ρ_parallel(ρ_A, ρ_O)
    ρ_A = ρ_A[1:end-1] .* ρ_A[2:end]
    ρ_O = ρ_O[1:end-1] .* ρ_O[2:end]
    ρ_A[1:2:end] .= NaN
    ρ_O[2:2:end] .= NaN
    return ρ_A, ρ_O
end
