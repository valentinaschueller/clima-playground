using Interpolations
using Plots
using LaTeXStrings
using clima_playground
import ClimaTimeSteppers as CTS
import ClimaCore as CC
import ClimaCoupler: Interfacer

function from_semtner_to_SI(flux)
    """
    Incoming unit: kcal / (cm^2 month). Outgoing: J / (m^2 s)
    """
    conversion_factor = 4184 * 1e4 / (3600 * 24 * 30)
    return flux .* conversion_factor
end

function interpolations()
    F_r = -[0, 0, 1.9, 9.9, 17.7, 19.2, 13.6, 9.0, 3.7, 0.4, 0, 0]
    F_L = -[10.4, 10.3, 10.3, 11.6, 15.1, 18.0, 19.1, 18.7, 16.5, 13.9, 11.2, 10.9]
    J_s = -[1.18, 0.76, 0.72, 0.29, -0.45, -0.39, -0.30, -0.40, -0.17, 0.10, 0.56, 0.79]
    J_q = [0.0, 0.02, 0.03, 0.09, 0.46, 0.70, 0.64, 0.66, 0.39, 0.19, 0.01, 0.01]

    F_r = from_semtner_to_SI(F_r)
    F_L = from_semtner_to_SI(F_L)
    J_s = from_semtner_to_SI(J_s)
    J_q = from_semtner_to_SI(J_q)

    months = 1:12
    F_r_spline = cubic_spline_interpolation(months, F_r; bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))
    F_L_spline = cubic_spline_interpolation(months, F_L; bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))
    J_s_spline = cubic_spline_interpolation(months, J_s; bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))
    J_q_spline = cubic_spline_interpolation(months, J_q; bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))

    seconds_in_year = range(0, 12, length=360)
    plot(seconds_in_year, F_r_spline(seconds_in_year); label=L"F_r")
    plot!(seconds_in_year, F_L_spline(seconds_in_year); label=L"F_L")
    plot!(seconds_in_year, J_s_spline(seconds_in_year); label=L"J_s")
    plot!(seconds_in_year, J_q_spline(seconds_in_year); label=L"J_q")
end

function ice_only_test()
    context = CC.ClimaComms.context()
    point_space = CC.Spaces.PointSpace(context, CC.Geometry.ZPoint(0.0))
    p = SimulationParameters(ice_model_type=:thickness_feedback, T_O_ini=273.15 - 1.9, Δu_IO=0.03)
    p.T_A = p.T_A_ini
    p.T_O = p.T_O_ini
    if p.ice_model_type != :constant
        @info("Determine initial ice surface temperature from SEB.")
        p.T_I_ini = solve_surface_energy_balance(p)[1]
    end
    odesolver = CTS.ExplicitAlgorithm(CTS.RK4())
    field_h_I = CC.Fields.ones(point_space) .* p.h_I_ini
    h_ice_0 = CC.Fields.FieldVector(data=field_h_I)
    output_dir = "ice_only"
    # rm(output_dir, recursive=true)
    mkpath(output_dir)
    ice_sim = ice_init(odesolver, h_ice_0, point_space, p, output_dir)
    Interfacer.step!(ice_sim, p.t_max)
end
