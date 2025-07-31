using Interpolations
using Plots
using LaTeXStrings
using Statistics
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

function interpolants()
    F_r = [0, 0, 1.9, 9.9, 17.7, 19.2, 13.6, 9.0, 3.7, 0.4, 0, 0]
    F_L = [10.4, 10.3, 10.3, 11.6, 15.1, 18.0, 19.1, 18.7, 16.5, 13.9, 11.2, 10.9]
    J_s = [1.18, 0.76, 0.72, 0.29, -0.45, -0.39, -0.30, -0.40, -0.17, 0.10, 0.56, 0.79]
    J_q = -[0.0, 0.02, 0.03, 0.09, 0.46, 0.70, 0.64, 0.66, 0.39, 0.19, 0.01, 0.01]

    F_r = from_semtner_to_SI(F_r)
    F_L = from_semtner_to_SI(F_L)
    J_s = from_semtner_to_SI(J_s)
    J_q = from_semtner_to_SI(J_q)

    months = 1:12
    F_r_spline = cubic_spline_interpolation(months, F_r; bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))
    F_L_spline = cubic_spline_interpolation(months, F_L; bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))
    J_s_spline = cubic_spline_interpolation(months, J_s; bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))
    J_q_spline = cubic_spline_interpolation(months, J_q; bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))

    return F_r_spline, F_L_spline, J_s_spline, J_q_spline
end

function plot_interpolants()
    F_r, F_L, J_s, J_q = interpolants()
    seconds_in_year = range(0, 12, length=360)
    plot(seconds_in_year, F_r(seconds_in_year); label=L"F_r")
    plot!(seconds_in_year, F_L(seconds_in_year); label=L"F_L")
    plot!(seconds_in_year, J_s(seconds_in_year); label=L"J_s")
    plot!(seconds_in_year, J_q(seconds_in_year); label=L"J_q")
end

function ice_only_test()
    conversion_factor = 4.184 * 1e7 / (3600 * 24 * 30 * 12)
    @info(conversion_factor)
    F_O = t -> 0.0 * conversion_factor
    context = CC.ClimaComms.context()
    point_space = CC.Spaces.PointSpace(context, CC.Geometry.ZPoint(0.0))
    F_r, F_L, J_s, J_q = interpolants()
    p = SimulationParameters(
        SW_in=F_r,
        LW_in=F_L,
        J_s=J_s,
        J_q=J_q,
        F_O=F_O,
        Δt_min=3600 * 24,
        h_I_ini=3.12,
        timestepping=:explicit,
        t_max=3600 * 24 * 30 * 12,
        ice_model_type=:thickness_feedback,
    )
    if p.ice_model_type != :constant
        @info("Determine initial ice surface temperature from SEB.")
        T_eq = compute_T_Is(p)
        p.T_I_ini = min(T_eq, 273)
    end
    odesolver = get_odesolver(Val(p.timestepping))
    field_h_I = CC.Fields.ones(point_space) .* p.h_I_ini
    h_ice_0 = CC.Fields.FieldVector(data=field_h_I)
    output_dir = "ice_only"
    rm(output_dir, recursive=true, force=true)
    mkpath(output_dir)
    ice_sim = ice_init(odesolver, h_ice_0, point_space, p, output_dir)
    Interfacer.step!(ice_sim, p.t_max)
    h = [parent(u)[1] for u in ice_sim.integrator.sol.u]
    t = ice_sim.integrator.sol.t ./ (3600 * 24 * 30)
    months = "JFMAMJJASONDJ"
    xticks = (0:12, months)
    plot(t, h, xticks=xticks)
    display(current())
    T_Is = [min(compute_T_Is(p, h[i], t[i]), 273) for i in range(1, length(t))]
    plot(t, T_Is .- 273, yflip=true, xticks=xticks)
    display(current())
    @info("Mean thickness: $(mean(h))")
    return ice_sim
end
