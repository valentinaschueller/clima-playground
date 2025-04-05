using Plots

function plot_obukhov_C_dependencies(C_variable_name)
    if !(C_variable_name in ["C_AO", "C_AI"])
        return nothing
    end
    physical_values = define_realistic_vals()
    a_is = C_variable_name == "C_AO" ? [physical_values[:a_i]] : 0:0.01:1

    pos_vars = -200:-50# L_AO or L_AI unstable
    neg_vars = 10:200# L_AO or L_AI stable
    C_pos = zeros(length(a_is), length(pos_vars))
    C_neg = zeros(length(a_is), length(neg_vars))
    for (i, a_i) in enumerate(a_is)
        for (j, pos_var) in enumerate(pos_vars)
            physical_values[:L_AO] =
                (C_variable_name == "C_AO") ? pos_var : physical_values[:L_AO]
            physical_values[:L_AI] =
                (C_variable_name == "C_AI") ? pos_var : physical_values[:L_AI]
            physical_values =
                (C_variable_name == "C_AO") ? update_C_AO(physical_values) :
                update_physical_values(a_i, physical_values)
            C_pos[i, j] =
                (C_variable_name == "C_AO") ? physical_values[:C_AO] :
                physical_values[:C_AI]
        end
        for (j, neg_var) in enumerate(neg_vars)
            physical_values[:L_AO] =
                (C_variable_name == "C_AO") ? neg_var : physical_values[:L_AO]
            physical_values[:L_AI] =
                (C_variable_name == "C_AI") ? neg_var : physical_values[:L_AI]
            physical_values =
                (C_variable_name == "C_AO") ? update_C_AO(physical_values) :
                update_physical_values(a_i, physical_values)
            C_neg[i, j] =
                (C_variable_name == "C_AO") ? physical_values[:C_AO] :
                physical_values[:C_AI]
        end
    end

    if C_variable_name == "C_AO"
        gr()
        println(maximum([maximum(C_pos), maximum(C_neg)]))
        println(minimum([minimum(C_pos), minimum(C_neg)]))
        plot(
            pos_vars,
            vec(C_pos),
            xlabel = L"$L^A_O$",
            ylabel = L"$C^A_O$",
            label = "",
            color = :black,
        )
        plot!(neg_vars, vec(C_neg), label = "", color = :black)
        display(current())
    elseif C_variable_name == "C_AI"
        plotly()
        println(maximum([maximum(C_pos), maximum(C_neg)]))
        println(minimum([minimum(C_pos), minimum(C_neg)]))
        surface(a_is, pos_vars, C_pos', xlabel = "aᴵ", ylabel = "Lᴬᴵ", zlabel = "Cᴬᴵ")
        surface!(a_is, neg_vars, C_neg')
    end
end

function analytical_convergence_factor_dependence()
    nus = range(0, stop = 10, length = 50)
    omegas = range(0.001, stop = 10, length = 50)
    a_is = range(0.001, stop = 1, length = 10)
    rhos = zeros(length(nus), length(omegas), length(a_is))
    physical_values = define_realistic_vals()

    # Loop over values
    for (k, a_i) in enumerate(a_is)
        for (i, nu) in enumerate(nus)
            for (j, omega) in enumerate(omegas)
                physical_values[:omega] = omega
                physical_values[:nu] = nu
                physical_values[:a_i] = a_i
                rho = compute_rho_with_nu_and_omega_variable(physical_values)
                rhos[i, j, k] = rho
            end
        end
        i_max, j_max = Tuple(CartesianIndices(rhos)[argmax(rhos)])
        nu_max = nus[i_max]
        omega_max = omegas[j_max]
        println("supremum at ν=$nu_max and ω=$omega_max")
    end


    plotly()
    surface(omegas, nus, rhos[:, :, 1], color = :viridis)
end

function compute_rho_with_nu_and_omega_variable(params)
    real_part_o = sqrt(
        (sqrt(params[:nu]^2 + params[:omega]^2) + params[:nu]) / (2 * params[:alpha_o]),
    )
    imag_part_o =
        im *
        sign(params[:omega]) *
        sqrt(
            (sqrt(params[:nu]^2 + params[:omega]^2) - params[:nu]) / (2 * params[:alpha_o]),
        )
    sigma_o = real_part_o + imag_part_o
    real_part_a = sqrt(
        (sqrt(params[:nu]^2 + params[:omega]^2) + params[:nu]) / (2 * params[:alpha_a]),
    )
    imag_part_a =
        im *
        sign(params[:omega]) *
        sqrt(
            (sqrt(params[:nu]^2 + params[:omega]^2) - params[:nu]) / (2 * params[:alpha_a]),
        )
    sigma_a = real_part_a + imag_part_a
    eta_AO =
        params[:C_AO] *
        abs(params[:u_atm] - params[:u_oce]) *
        params[:rho_atm] *
        params[:c_atm]
    eta_OI = params[:C_OI] * abs(params[:u_oce]) * params[:rho_oce] * params[:c_oce]
    eta_AI = params[:C_AI] * abs(params[:u_atm]) * params[:rho_atm] * params[:c_atm]
    return abs(
        (1 - params[:a_i])^2 * eta_AO^2 / (
            (
                params[:k_oce] *
                sigma_o *
                (1 / tanh(sigma_o * (params[:h_oce] - params[:z_0numO]))) +
                (1 - params[:a_i]) * eta_AO +
                params[:a_i] * eta_OI
            ) * (
                params[:k_atm] *
                sigma_a *
                (1 / tanh(sigma_a * (params[:h_atm] - params[:z_0numA]))) +
                (1 - params[:a_i]) * eta_AO +
                params[:a_i] * eta_AI
            )
        ),
    )
end
