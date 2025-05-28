import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities

export get_vertical_space, get_coupled_sim

function get_vertical_space(device, lower_boundary, upper_boundary, nelems)
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(lower_boundary),
        CC.Geometry.ZPoint{Float64}(upper_boundary);
        boundary_names=(:bottom, :top),
    )
    mesh = CC.Meshes.IntervalMesh(domain, nelems=nelems)
    return CC.Spaces.CenterFiniteDifferenceSpace(device, mesh)
end

function get_stable_range(initial_value_fvs)
    max_value = floatmin()
    min_value = floatmax()
    for fv in initial_value_fvs
        min_value = min(min_value, minimum(fv.data))
        max_value = max(max_value, maximum(fv.data))
    end
    return min_value - eps(min_value), max_value + eps(max_value)
end

function get_cit_boundary_space(device, p::SimulationParameters)
    time_points = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint(p.t_0),
        CC.Geometry.ZPoint(p.Δt_cpl);
        boundary_names=(:start, :end),
    )
    time_mesh = CC.Meshes.IntervalMesh(
        time_points,
        nelems=Int(p.timerange[end] / p.Δt_min),
    )
    return CC.Spaces.FaceFiniteDifferenceSpace(device, time_mesh)
end

function get_coupled_sim(p::SimulationParameters)
    context = CC.ClimaComms.context()
    device = CC.ClimaComms.device(context)
    output_dir = "output"
    rm(output_dir, recursive=true)
    mkpath(output_dir)
    dir_paths = (
        output=output_dir,
        artifacts=output_dir,
        regrid=output_dir,
        checkpoints=output_dir,
    )

    center_space_atm = get_vertical_space(device, p.z_A0, p.h_A, p.n_A)
    center_space_oce = get_vertical_space(device, -p.h_O, -p.z_O0, p.n_O)
    point_space = CC.Spaces.PointSpace(context, CC.Geometry.ZPoint(0.0))

    field_atm = CC.Fields.ones(center_space_atm) .* p.T_A_ini
    field_oce = CC.Fields.ones(center_space_oce) .* p.T_O_ini
    T_atm_0 = CC.Fields.FieldVector(data=field_atm)
    T_oce_0 = CC.Fields.FieldVector(data=field_oce)

    if p.ice_model_type != :constant
        @info("Determine initial ice surface temperature from SEB.")
        cache = (; p_dict..., T_A=T_atm_0[1] .* CC.Fields.ones(point_space))
        p.T_I_ini = solve_surface_energy_balance(cache)[1]
    end
    field_ice = CC.Fields.ones(point_space) .* p.T_I_ini
    field_h_I = CC.Fields.ones(point_space) .* p.h_I_ini
    T_ice_0 = CC.Fields.FieldVector(data=field_ice)
    h_ice_0 = CC.Fields.FieldVector(data=field_h_I)

    p.stable_range = get_stable_range([T_atm_0, T_oce_0, (data=[p.T_I_ini, p.T_Ib],)])
    if p.ice_model_type == :thickness_feedback
        p.stable_range = nothing
    end

    if p.boundary_mapping == "cit"
        boundary_space = get_cit_boundary_space(device, p)
    else
        boundary_space = point_space
    end
    p.T_O = T_oce_0[end] .* CC.Fields.ones(boundary_space)
    p.T_A = T_atm_0[1] .* CC.Fields.ones(boundary_space)
    p.T_Is = T_ice_0[1] .* CC.Fields.ones(boundary_space)

    odesolver = CTS.ExplicitAlgorithm(CTS.RK4())
    atmos_sim = atmos_init(odesolver, T_atm_0, center_space_atm, p, output_dir)
    ocean_sim = ocean_init(odesolver, T_oce_0, center_space_oce, p, output_dir)
    ice_sim = ice_init(odesolver, h_ice_0, point_space, p, output_dir)

    start_date = "19790301"
    date = Dates.DateTime(start_date, Dates.dateformat"yyyymmdd")
    dates = (;
        date=[date],
        date0=[date],
        date1=[Dates.firstdayofmonth(date)],
        new_month=[false],
    )

    model_sims = (atmos_sim=atmos_sim, ocean_sim=ocean_sim, ice_sim=ice_sim)

    coupler_field_names = []
    for sim in model_sims
        Interfacer.add_coupler_fields!(coupler_field_names, sim)
    end
    coupler_fields = Interfacer.init_coupler_fields(Float64, coupler_field_names, boundary_space)

    cs = Interfacer.CoupledSimulation{Float64}(
        context,
        dates,
        boundary_space,
        coupler_fields,
        nothing, # conservation checks
        (p.t_0, p.t_max),
        p.Δt_cpl,
        model_sims,
        (;), # callbacks
        dir_paths,
        FluxCalculator.PartitionedStateFluxes(),
        nothing, # thermo_params
        nothing, # diagnostic_handlers
    )
    return cs
end
