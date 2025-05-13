export is_stable, has_converged, extract_ρ, update_ρ_parallel, UnstableError

struct UnstableError <: Exception
end

function is_stable(values, value_range)
    if (
        any(isnan, values) ||
        maximum(values) > maximum(value_range) ||
        minimum(values) < minimum(value_range)
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
