"""
Extracts an Array from an Array of FieldVectors.

**Arguments:**

-`field_vecs: Array`: An array of fieldvectors containing temperature values.
-`domain: String`: Either `"atm"` or `"oce"`.

"""
function extract_matrix(field_vecs, domain)
    matrix = []
    for field_vec in field_vecs
        if domain == "atm"
            field = field_vec.atm
            values = parent(field)
            push!(matrix, values)
        else
            field = field_vec.oce
            values = parent(field)
            push!(matrix, values)
        end
    end
    return hcat(matrix...)
end

"""
Checks if the computed temperatures are reasonable or if the model has gone unstable.

**Arguments:**

-`atmos_vals: Array`: Atmosphere temperatures.
-`ocean_vals: Array`: Ocean temperatures.
-`upper_limit_temp: Float64`: Maximum reasonable temperature.
-`lower_limit_temp: Float64`: Minimum reasonable temperature.
-`iter: Int`: Current iteration.

"""
function is_stable(atmos_vals, ocean_vals, upper_limit_temp, lower_limit_temp, iter)
    stable = true
    stopped_at_nan_atm = false
    stopped_at_nan_oce = false
    if (
        any(isnan, ocean_vals) ||
        maximum(ocean_vals) > upper_limit_temp ||
        minimum(ocean_vals) < lower_limit_temp
    )
        println("stopped due to instability in ocean model")
        stable = false
        if iter == 1
            stopped_at_nan_oce = true
        end
    elseif (
        any(isnan, atmos_vals) ||
        maximum(atmos_vals) > upper_limit_temp ||
        minimum(atmos_vals) < lower_limit_temp
    )
        println("stopped due to instability in atmosphere model")
        stable = false
        if iter == 1
            stopped_at_nan_atm = true
        end
    end
    return stable, stopped_at_nan_atm, stopped_at_nan_oce
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
    pre_bound_ocean_vals,
    iter,
)
    if iter == 1
        return false
    end
    bound_errors_atm_iter = abs.(bound_atmos_vals .- pre_bound_atmos_vals)
    bound_errors_oce_iter = abs.(bound_ocean_vals .- pre_bound_ocean_vals)
    tols_atm = 100 * eps.(max.(abs.(bound_atmos_vals), abs.(pre_bound_atmos_vals)))
    tols_oce = 100 * eps.(max.(abs.(bound_ocean_vals), abs.(pre_bound_ocean_vals)))
    if all(bound_errors_atm_iter .< tols_atm)
        println("Stopped at iter $iter for the atmosphere")
        return true
    elseif all(bound_errors_oce_iter .< tols_oce)
        println("Stopped at iter $iter for the ocean")
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
-`conv_facs_oce::Array`: Ocean convergence factors.
-`conv_facs_atm::Array`: Atmosphere convergence factors.
-`physical_values::Dict`: Can be defined using `define_realistic_vals()`.

**Optional Keyword Arguments:**

-`dims::Int`: Dimensions over which to  reverse the convergence factor matrices, default: `1`.
-`param_analytic::Array`: The analytical values for the variable, default: `nothing`.
-`conv_facs_analytic::Array`: The analytical convergence factors, default: `nothing`.

"""
function handle_variable(
    var,
    var_name,
    conv_facs_oce,
    conv_facs_atm,
    physical_values;
    dims=1,
    param_analytic=nothing,
    conv_facs_analytic=nothing,
)
    # Special treatment of n_atm and n_oce.
    if var_name == "n_atm"
        var = (physical_values[:h_atm] - physical_values[:z_0numA]) ./ reverse(var)
        if !isnothing(conv_facs_analytic) && !isnothing(param_analytic)
            param_analytic = reverse(
                (physical_values[:h_atm] - physical_values[:z_0numA]) ./ param_analytic,
            )
        end
    elseif var_name == "n_oce"
        var = (physical_values[:h_oce] - physical_values[:z_0numO]) ./ reverse(var)
        if !isnothing(conv_facs_analytic) && !isnothing(param_analytic)
            param_analytic = reverse(
                (physical_values[:h_oce] - physical_values[:z_0numO]) ./ param_analytic,
            )
        end
    elseif var_name == "n_t_atm" || var_name == "n_t_oce"
        var = physical_values[:Δt_min] ./ reverse(var)
        if !isnothing(conv_facs_analytic) && !isnothing(param_analytic)
            param_analytic = reverse(physical_values[:Δt_min] ./ param_analytic)
        end
    end
    if var_name in ["n_atm", "n_oce", "n_t_atm", "n_t_oce"]
        conv_facs_oce =
            !isnothing(conv_facs_oce) ? reverse(conv_facs_oce, dims=dims) : nothing
        conv_facs_atm =
            !isnothing(conv_facs_atm) ? reverse(conv_facs_atm, dims=dims) : nothing
        conv_facs_analytic =
            !isnothing(conv_facs_analytic) ? reverse(conv_facs_analytic, dims=dims) :
            nothing
    end
    return var, conv_facs_oce, conv_facs_atm, param_analytic, conv_facs_analytic
end

""" 
Extracts the first convergence factors.

**Arguments:**

-`conv_fac_atm::Array`: Atmosphere convergence factor.
-`conv_fac_oce::Array`: Ocean convergence factor.

"""
function extract_conv_fac(conv_fac_atm, conv_fac_oce)
    if conv_fac_atm isa AbstractArray
        conv_fac_atm = !isempty(conv_fac_atm) ? conv_fac_atm[1] : NaN
    end
    if conv_fac_oce isa AbstractArray
        conv_fac_oce = !isempty(conv_fac_oce) ? conv_fac_oce[1] : NaN
    end
    return conv_fac_atm, conv_fac_oce
end
