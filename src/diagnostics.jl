export get_prognostic_data!

function get_prognostic_data!(out, Y, p, t)
    if isnothing(out)
        return copy(Y.data)
    else
        out .= Y.data
    end
end

function get_diagnostic(variable, space, Δt, output_dir)
    output_schedule = CD.Schedules.EveryDtSchedule(Δt)
    netcdf_writer = CD.Writers.NetCDFWriter(space, output_dir)
    scheduled_temperature = CD.ScheduledDiagnostic(
        variable=variable,
        output_writer=netcdf_writer,
        compute_schedule_func=output_schedule,
    )

    return scheduled_temperature
end
