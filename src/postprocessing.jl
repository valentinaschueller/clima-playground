export check_stability, has_converged, compute_ϱ_numerical, UnstableError

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
            @warn "Outside initial value range!"
        end
    end
end

function has_converged(T_A_Γ, T_A_Γ_old, T_O_Γ, T_O_Γ_old)
    if isnothing(T_A_Γ_old)
        return false
    end
    e_A_Γ = abs.(T_A_Γ .- T_A_Γ_old)
    e_O_Γ = abs.(T_O_Γ .- T_O_Γ_old)
    tols_atm = 100 * eps.(max.(abs.(T_A_Γ), abs.(T_A_Γ_old)))
    tols_oce = 100 * eps.(max.(abs.(T_O_Γ), abs.(T_O_Γ_old)))
    if all(e_A_Γ .< tols_atm) && all(e_O_Γ .< tols_oce)
        return true
    else
        return false
    end
end

function compute_ϱ_numerical(coupling_variable, parallel::Bool)
    ϱ = []
    e_old = abs.(coupling_variable[1] .- coupling_variable[end])
    for i = 2:length(coupling_variable)-1
        e_new = abs.(coupling_variable[i] .- coupling_variable[end])
        ϱ_i = norm(e_new) / norm(e_old)
        push!(ϱ, ϱ_i)

        e_old = e_new
    end
    if parallel
        update_ϱ_parallel!(ϱ)
    end
    return mean_ϱ(ϱ)
end

function mean_ϱ(ϱ_array)
    non_nan_values = filter(!isnan, ϱ_array)
    if isempty(non_nan_values)
        return NaN
    end
    return mean(non_nan_values)
end

function update_ϱ_parallel!(ϱ)
    ϱ = ϱ[1:end-1] .* ϱ[2:end]
    ϱ[1:2:end] .= NaN
end
