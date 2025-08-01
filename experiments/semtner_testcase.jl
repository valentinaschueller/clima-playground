using Interpolations
using Plots
using LaTeXStrings
using Statistics
using NetCDF
using clima_playground
import ClimaTimeSteppers as CTS
import ClimaCore as CC
import ClimaCoupler: Interfacer
import ClimaDiagnostics as CD

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

    halfmonth = 15 * 24 * 3600
    month = 2 * halfmonth
    jan15 = halfmonth
    dec15 = 11 * month + halfmonth
    seconds = jan15:month:dec15
    F_r_spline = cubic_spline_interpolation(seconds, F_r; bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))
    F_L_spline = cubic_spline_interpolation(seconds, F_L; bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))
    J_s_spline = cubic_spline_interpolation(seconds, J_s; bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))
    J_q_spline = cubic_spline_interpolation(seconds, J_q; bc=Periodic(OnCell()), extrapolation_bc=Periodic(OnCell()))

    return F_r_spline, F_L_spline, J_s_spline, J_q_spline
end

function plot_interpolants()
    F_r, F_L, J_s, J_q = interpolants()
    time_in_s = range(0, 24 * 3600 * 30 * 12, length=360)
    time_in_months = time_in_s ./ (24 * 3600 * 30)
    plot(time_in_months, F_r(time_in_s); label=L"F_r")
    plot!(time_in_months, F_L(time_in_s); label=L"F_L")
    plot!(time_in_months, J_s(time_in_s); label=L"J_s")
    plot!(time_in_months, J_q(time_in_s); label=L"J_q")
end

function ice_only_test()
    conversion_factor = 4.184 * 1e7 / (3600 * 24 * 30 * 12)
    @info(conversion_factor)
    F_O = t -> 1.5 * conversion_factor
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
        h_I_ini=2.65,
        timestepping=:explicit,
        t_max=3600 * 24 * 30 * 12,
        ice_model_type=:thickness_feedback,
        alb_I=0.8,
    )
    if p.ice_model_type != :constant
        @info("Determine initial ice surface temperature from SEB.")
        p.T_I_ini = T_Is(p)
    end
    odesolver = get_odesolver(Val(p.timestepping))
    field_h_I = CC.Fields.ones(point_space) .* p.h_I_ini
    h_ice_0 = CC.Fields.FieldVector(data=field_h_I)
    output_dir = "ice_only"
    rm(output_dir, recursive=true, force=true)
    mkpath(output_dir)
    ice_sim = ice_init(odesolver, h_ice_0, point_space, p, output_dir)
    Interfacer.step!(ice_sim, p.t_max)
    dt = CD.seconds_to_str_short(p.Δt_min)
    time = ncread("$(output_dir)/h_I_$(dt)_inst.nc", "time")
    h_I = ncread("$(output_dir)/h_I_$(dt)_inst.nc", "h_I")
    T_Is = ncread("$(output_dir)/T_Is_$(dt)_inst.nc", "T_Is") .- 273
    t = time ./ (3600 * 24 * 30)
    months = "JFMAMJJASONDJ"
    xticks = (0:12, months)
    p1 = plot(t, h_I, xticks=xticks, label=L"h_I", ylabel="Ice Thickness [m]", color=:black)
    p2 = plot(t, T_Is, yflip=true, xticks=xticks, label=L"T_{I,s}", ylabel="Temperature [°C]", xlabel="Month", color=:black, legend=:bottomright)
    plot(p1, p2, layout=(2, 1), legendfontsize=12, linewidth=2)
    display(current())
    @info("Mean thickness: $(mean(h_I))")
    return ice_sim
end
