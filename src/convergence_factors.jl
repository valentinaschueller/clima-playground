import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer,
    FieldExchanger,
    FluxCalculator,
    Interfacer,
    TimeManager,
    Utilities

function compute_actual_rho(physical_values)
    sigma_o = sqrt(abs(physical_values[:w_min]) / (2 * physical_values[:alpha_o])) * (1 + im)
    sigma_a = sqrt(abs(physical_values[:w_min]) / (2 * physical_values[:alpha_a])) * (1 + im)
    eta_AO = physical_values[:C_AO] * abs(physical_values[:u_atm] - physical_values[:u_oce]) * physical_values[:rho_atm] * physical_values[:c_atm]
    eta_OI = physical_values[:C_OI] * abs(physical_values[:u_oce]) * physical_values[:rho_oce] * physical_values[:c_oce]
    eta_AI = physical_values[:C_AI] * abs(physical_values[:u_atm]) * physical_values[:rho_atm] * physical_values[:c_atm]
    return abs((1 - physical_values[:a_i])^2 * eta_AO^2 / ((physical_values[:k_oce] * sigma_o * (1 / tanh(sigma_o * (physical_values[:h_oce] - physical_values[:z_0numO]))) + (1 - physical_values[:a_i]) * eta_AO + physical_values[:a_i] * eta_OI) * (physical_values[:k_atm] * sigma_a * (1 / tanh(sigma_a * (physical_values[:h_atm] - physical_values[:z_0numA]))) + (1 - physical_values[:a_i]) * eta_AO + physical_values[:a_i] * eta_AI)))
end

function compute_numerical_conv_fac(physical_values; return_conv_facs=true, plot_conv_facs=false, print_conv_facs=false)
    cs = get_coupled_sim(physical_values)
    if return_conv_facs
        conv_fac_atm, conv_fac_oce = solve_coupler!(cs, print_conv=false, plot_conv=false, return_conv=return_conv_facs)
        return conv_fac_atm, conv_fac_oce
    elseif plot_conv_facs || print_conv_facs
        solve_coupler!(cs, print_conv=print_conv_facs, plot_conv=plot_conv_facs, return_conv=return_conv_facs)
    end

end

function get_conv_facs_one_variable(physical_values, var2s, var2_name; a_i_variable=nothing, analytic=false, log_scale=false)
    # Introducing the convergence factors
    vary_a_i = !isnothing(a_i_variable)
    a_i_variable = vary_a_i ? a_i_variable : [physical_values[:a_i]]
    conv_facs_atm = zeros(length(a_i_variable), length(var2s))
    conv_facs_oce = zeros(length(a_i_variable), length(var2s))

    # Create range if we want to calculate the analytical convergence factor
    if analytic
        if log_scale
            variable2_range = exp10.(range(log10(var2s[1]), log10(var2s[end]), length=100))
        else
            variable2_range = range(var2s[1], var2s[end], length=100)
        end
        conv_facs_analytic = zeros(length(a_i_variable), length(variable2_range))
    else
        variable2_range = nothing
    end

    if var2_name == "h_atm"
        delta_z_atm = (physical_values[:h_atm] - physical_values[:z_0numA]) / physical_values[:n_atm]
    elseif var2_name == "h_oce"
        delta_z_oce = (physical_values[:h_oce] - physical_values[:z_0numO]) / physical_values[:n_oce]
    end

    # Looping over the variables
    for (j, a_i) in enumerate(a_i_variable)
        physical_values[:a_i] = a_i
        physical_values = update_physical_values(a_i, physical_values)
        for (k, var2) in enumerate(var2s)
            # Update physical values based on variable 1
            physical_values[Symbol(var2_name)] = var2
            if var2_name == "a_i"
                physical_values = update_physical_values(var2, physical_values)
            elseif var2_name == "h_atm"
                physical_values[:n_atm] = Int((var2 - physical_values[:z_0numA]) / delta_z_atm)
            elseif var2_name == "h_oce"
                physical_values[:n_oce] = Int((var2 - physical_values[:z_0numO]) / delta_z_oce)
            elseif var2_name == "t_max"
                physical_values[:delta_t_cpl] = var2
            end
            conv_fac_atm, conv_fac_oce = compute_numerical_conv_fac(physical_values)
            conv_facs_atm[j, k], conv_facs_oce[j, k] = extract_conv_fac(conv_fac_atm, conv_fac_oce)
        end
        if analytic
            for (k, var2) in enumerate(variable2_range)
                physical_values[Symbol(var2_name)] = var2
                if var2_name == "a_i"
                    physical_values = update_physical_values(var2, physical_values)
                elseif var2_name == "t_max"
                    physical_values[:w_min] = pi / var2
                end
                conv_fac_analytic = compute_actual_rho(physical_values)
                conv_facs_analytic[j, k] = conv_fac_analytic
            end
        else
            conv_facs_analytic = nothing
        end
    end
    return conv_facs_atm, conv_facs_oce, variable2_range, conv_facs_analytic
end

function stability_check(physical_values, n_zs, n_ts, var1_name, var2_name)
    domain = var1_name[end-2:end]
    unstable_matrix = zeros(length(n_ts), length(n_zs))
    theoretical_vals_matrix = fill((0.0, 0.0, 0.0), length(n_ts), length(n_zs))

    for (i, n_z) in enumerate(n_zs)
        for (j, n_t) in enumerate(n_ts)
            physical_values[Symbol(var1_name)] = n_z
            physical_values[Symbol(var2_name)] = n_t

            cs = get_coupled_sim(physical_values)

            starting_temp_oce = cs.model_sims.atmos_sim.params.T_atm_ini
            starting_temp_atm = cs.model_sims.ocean_sim.params.T_oce_ini
            starting_temp_ice = cs.model_sims.atmos_sim.params.T_ice_ini
            upper_limit_temp = maximum([starting_temp_oce, starting_temp_atm, starting_temp_ice])
            lower_limit_temp = minimum([starting_temp_oce, starting_temp_atm, starting_temp_ice])

            delta_t = physical_values[:delta_t_min] ./ n_t

            if domain == "atm"
                Interfacer.step!(cs.model_sims.atmos_sim, physical_values[:delta_t_cpl])
                states = copy(cs.model_sims.atmos_sim.integrator.sol.u)

                delta_z = (physical_values[:h_atm] - physical_values[:z_0numA]) ./ n_z
                theoretical_vals_matrix[i, j] = (2 * physical_values[:k_oce] * delta_t / (physical_values[:c_oce] * physical_values[:rho_oce] * delta_z^2),
                    2 * (1 - physical_values[:a_i]) * physical_values[:C_AO] * abs(physical_values[:u_atm] - physical_values[:u_oce]) * delta_t ./ delta_z,
                    2 * physical_values[:a_i] * physical_values[:C_AI] * abs(physical_values[:u_atm]) * delta_t ./ delta_z
                )
            elseif domain == "oce"
                Interfacer.step!(cs.model_sims.ocean_sim, physical_values[:delta_t_cpl])
                states = copy(cs.model_sims.ocean_sim.integrator.sol.u)

                # TODO: This is not correct
                delta_z = (physical_values[:h_oce] - physical_values[:z_0numO]) ./ n_z
                theoretical_vals_matrix[i, j] = (2 * physical_values[:k_oce] * delta_t / (physical_values[:c_oce] * physical_values[:rho_oce] * delta_z^2),
                    2 * (1 - physical_values[:a_i]) * physical_values[:C_AO] * abs(physical_values[:u_atm] - physical_values[:u_oce]) * delta_t ./ delta_z,
                    2 * physical_values[:a_i] * physical_values[:C_AI] * abs(physical_values[:u_atm]) * delta_t ./ delta_z
                )
            end

            vals = extract_matrix(states, domain)
            if (any(isnan, vals) || maximum(vals) > upper_limit_temp || minimum(vals) < lower_limit_temp)
                println("unstable")
                unstable_matrix[i, j] = Inf
            else
                unstable_matrix[i, j] = NaN
            end
        end
    end
    return unstable_matrix, theoretical_vals_matrix
end