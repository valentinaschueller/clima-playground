import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities

function psi(xi, atm, stable, heat)
    if stable
        if atm && heat
            return -2 / 3 * (xi - (5 / 0.35)) * exp(-0.35 * xi) - (1 + (2 * xi / 3))^1.5 -
                   ((2 / 3) * 5 / 0.35) + 1
        elseif atm && !heat
            return -2 / 3 * (xi - (5 / 0.35)) * exp(-0.35 * xi) - xi - ((2 / 3) * 5 / 0.35)
        elseif !atm
            return -5 * xi
        end
    else
        if !atm && heat
            x = (1 - 25 * xi)^(1 / 3)
            return sqrt(3) * (atan(sqrt(3)) - atan(1 / sqrt(3) * (2 * x + 1))) +
                   (3 / 2) * log((x^2 + x + 1) / 3)
        elseif !atm && !heat
            x = (1 - 14 * xi)^(1 / 3)
            return sqrt(3) * (atan(sqrt(3)) - atan(1 / sqrt(3) * (2 * x + 1))) +
                   (3 / 2) * log((x^2 + x + 1) / 3)
        elseif atm && heat
            return 2 * log((1 + (1 - 16 * xi)^(1 / 2)) / 2)
        else
            x = (1 - 16 * xi)^(1 / 4)
            return pi / 2 - 2 * atan(x) + log((1 + x)^2 * (1 + x^2) / 8)
        end
    end
end

function define_realistic_vals()

    physical_values = Dict(
        :a_i => Float64(0),
        :rho_atm => Float64(1.225),
        :rho_oce => Float64(1000.0),
        :c_atm => Float64(1005.0),
        :c_oce => Float64(4182.0),
        :lambda_u => Float64(sqrt(1.225 / 1000.0)),
        :lambda_T => Float64(sqrt(1.225 / 1000.0) * 1005.0 / 4182.0),
        :nu_O => Float64(1e-6),
        :nu_A => Float64(1.5e-5),
        :mu => Float64(1e-6 / 1.5e-5),
        :kappa => Float64(0.4),
        :k_atm => Float64(0.02364),
        :k_oce => Float64(0.58),
        :alpha_o => Float64(0.58 / (1000.0 * 4182.0)),
        :alpha_a => Float64(0.02364 / (1.225 * 1005.0)),
        :alpha_eos => Float64(1.8e-4),
        :z_0numA => Float64(10.0),
        :z_0numO => Float64(1.0),
        :z_ruAO => Float64(2e-4),
        :z_ruAI => nothing,
        :z_rTAO => Float64(2e-4),
        :z_rTAI => Float64(1e-3),
        :C_OI => Float64(5e-3),
        :h_oce => Float64(51.0),
        :h_atm => Float64(210.0),
        :u_atm => Float64(5.0),
        :u_oce => Float64(1.0),
        :L_AO => Float64(50.0),
        :L_AI => Float64(50.0),
        :L_OA => nothing,
        :C_AO => nothing,
        :C_AI => nothing,
        :T_atm_ini => Float64(267.0),
        :T_oce_ini => Float64(271.0),
        :T_ice_ini => Float64(270.0),
        :t_max => Float64(3600.0),
        :delta_t_cpl => Float64(100),
        :delta_t_min => Float64(1.0),
        :n_t_atm => 50,
        :n_t_oce => 1,
        :n_atm => 200,
        :n_oce => 50,
        :boundary_mapping => "mean",
        :sin_field_atm => false,
        :sin_field_oce => false,
    )
    physical_values = update_C_AO(physical_values)
    physical_values[:w_min] = pi / physical_values[:t_max]
    z_ruAI = Float64(
        max(
            1e-3,
            0.93e-3 * (1 - physical_values[:a_i]) +
            6.05e-3 * exp(-17 * (physical_values[:a_i] - 0.5)^2),
        ),
    )
    physical_values[:z_ruAI] = z_ruAI
    physical_values[:C_AI] =
        physical_values[:kappa]^2 / (
            (
                log(physical_values[:z_0numA] / physical_values[:z_ruAI]) - psi(
                    physical_values[:z_0numA] / physical_values[:L_AI],
                    true,
                    physical_values[:L_AI] > 0,
                    false,
                )
            ) * (
                log(physical_values[:z_0numA] / physical_values[:z_rTAI]) - psi(
                    physical_values[:z_0numA] / physical_values[:L_AI],
                    true,
                    physical_values[:L_AI] > 0,
                    true,
                )
            )
        )
    return physical_values
end

function update_physical_values(a_i, physical_values)
    physical_values[:a_i] = a_i
    z_ruAI = Float64(max(1e-3, 0.93e-3 * (1 - a_i) + 6.05e-3 * exp(-17 * (a_i - 0.5)^2)))
    physical_values[:z_ruAI] = z_ruAI
    physical_values[:C_AI] =
        physical_values[:kappa]^2 / (
            (
                log(physical_values[:z_0numA] / physical_values[:z_ruAI]) - psi(
                    physical_values[:z_0numA] / physical_values[:L_AI],
                    true,
                    physical_values[:L_AI] > 0,
                    false,
                )
            ) * (
                log(physical_values[:z_0numA] / physical_values[:z_rTAI]) - psi(
                    physical_values[:z_0numA] / physical_values[:L_AI],
                    true,
                    physical_values[:L_AI] > 0,
                    true,
                )
            )
        )
    return physical_values
end

function update_C_AO(physical_values)
    physical_values[:L_OA] =
        physical_values[:lambda_u]^2 / (
            physical_values[:T_atm_ini] *
            physical_values[:alpha_eos] *
            physical_values[:lambda_T]
        )
    physical_values[:C_AO] =
        physical_values[:kappa]^2 / (
            (
                log(physical_values[:z_0numA] / physical_values[:z_ruAO]) - psi(
                    physical_values[:z_0numA] / physical_values[:L_AO],
                    true,
                    physical_values[:L_AO] > 0,
                    false,
                ) +
                physical_values[:lambda_u] * (
                    log(
                        physical_values[:lambda_u] * physical_values[:z_0numO] /
                        (physical_values[:z_ruAO] * physical_values[:mu]),
                    ) - psi(
                        physical_values[:z_0numO] / physical_values[:L_OA],
                        false,
                        physical_values[:L_OA] > 0,
                        false,
                    )
                )
            ) * (
                log(physical_values[:z_0numA] / physical_values[:z_rTAO]) - psi(
                    physical_values[:z_0numA] / physical_values[:L_AO],
                    true,
                    physical_values[:L_AO] > 0,
                    true,
                ) +
                physical_values[:lambda_T] * (
                    log(
                        physical_values[:lambda_T] * physical_values[:z_0numO] /
                        (physical_values[:z_rTAO] * physical_values[:mu]),
                    ) - psi(
                        physical_values[:z_0numO] / physical_values[:L_OA],
                        false,
                        physical_values[:L_OA] > 0,
                        true,
                    )
                )
            )
        )
    return physical_values
end

function get_var_dict()
    variable_dict = Dict(
        :t_max =>
            [Float64.([10, 100, 1000, 10000, 100000, 1000000]), L"$\Delta t_{cpl}$"],
        :h_atm => [Float64.([15, 50, 100, 200, 300, 500, 700, 1500]), L"$H^A$"],
        :h_oce => [Float64.([6, 10, 50, 100, 125, 150, 175, 200]), L"$H^O$"],
        :n_atm => [[20, 200, 2000, 20000, 200000], L"$\Delta z^A$"],
        :n_oce => [[5, 50, 500, 5000, 50000], L"$\Delta z^O$"],
        :C_AI => [
            Float64.(LinRange(0.0008691059360985882, 0.0027815836784748733, 10)),
            L"$C^A_I$",
        ],
        :C_AO => [
            Float64.(LinRange(0.0006765386884900067, 0.0014898726741724184, 10)),
            L"$C^A_O$",
        ],
        :C_OI => [Float64.(LinRange(0.001, 0.01, 10)), L"$C^O_I$"],
        :u_atm => [Float64.([0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]), L"$u^A$"],
        :u_oce => [Float64.([0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]), L"$u^O$"],
        :T_atm_ini =>
            [Float64.([260, 262, 264, 266, 268, 270, 273, 276]), L"$T^A(0,z)$"],
        :T_oce_ini => [Float64.([270, 271, 272, 273, 274, 275, 276]), L"$T^O(0,z)$"],
        :T_ice_ini => [Float64.([260, 262, 264, 266, 268, 270, 273]), L"$T^I$"],
        :n_t_atm =>
            [Float64.(10 .^ LinRange(log10(1), log10(1000), 10)), L"$\Delta t^A$"],
        :n_t_oce =>
            [Float64.(10 .^ LinRange(log10(1), log10(1000), 10)), L"$\Delta t^O$"],
    )
    return variable_dict
end

function get_color_dict()
    color_dict = Dict(
        0.0 => RGB(0.850, 0.325, 0.098),
        0.1 => RGB(0.0, 0.447, 0.741),
        0.2 => RGB(0.929, 0.694, 0.125),
        0.3 => RGB(0.494, 0.184, 0.556),
        0.4 => RGB(0.75, 0.0, 0.75),
        0.5 => RGB(0.301, 0.745, 0.933),
        0.6 => RGB(0.5, 0.5, 0.5),
        0.7 => RGB(0.466, 0.674, 0.188),
        0.8 => RGB(1.0, 0.843, 0.0),
        0.9 => RGB(0.0, 0.75, 0.75),
        1.0 => RGB(0.635, 0.078, 0.184),
    )
    linestyle_dict = Dict(
        0.0 => :solid,
        0.1 => :solid,
        0.2 => :dashdot,
        0.3 => :dash,
        0.4 => :dash,
        0.5 => :dashdotdot,
        0.6 => :solid,
        0.7 => :dot,
        0.8 => :dashdot,
        0.9 => :dashdotdot,
        1.0 => :dot,
    )
    return color_dict, linestyle_dict
end
