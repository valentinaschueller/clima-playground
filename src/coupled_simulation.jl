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

function initial_value_range(initial_value_fvs)
    max_value = floatmin()
    min_value = floatmax()
    for fv in initial_value_fvs
        min_value = min(min_value, minimum(fv.data))
        max_value = max(max_value, maximum(fv.data))
    end
    return min_value, max_value
end

function get_coupled_sim(p::SimulationParameters)
    context = CC.ClimaComms.context()
    device = CC.ClimaComms.device(context)

    center_space_atm = get_vertical_space(device, p.z_A0, p.h_A, p.n_A)
    center_space_oce = get_vertical_space(device, -p.h_O, -p.z_O0, p.n_O)

    coord = CC.Geometry.ZPoint{Float64}(0.0)
    point_space_ice = CC.Spaces.PointSpace(context, coord)

    atmos_facespace = CC.Spaces.FaceFiniteDifferenceSpace(center_space_atm)
    boundary_space = CC.Spaces.level(
        atmos_facespace,
        CC.Utilities.PlusHalf(CC.Spaces.nlevels(atmos_facespace) - 1),
    )

    stepping = (;
        Δt_min=Float64(p.Δt_min),
        timerange=(Float64(0.0), Float64(p.t_max)),
        Δt_coupler=Float64(p.Δt_cpl),
        odesolver=CTS.ExplicitAlgorithm(CTS.RK4()),
        nsteps_atm=p.n_t_A,
        nsteps_oce=p.n_t_O,
        nsteps_ice=p.n_t_I,
    )

    field_atm = CC.Fields.ones(Float64, center_space_atm) .* p.T_A_ini
    field_oce = CC.Fields.ones(Float64, center_space_oce) .* p.T_O_ini
    field_ice = CC.Fields.ones(Float64, point_space_ice) .* p.T_I_ini
    field_h_I = CC.Fields.ones(Float64, point_space_ice) .* p.h_I_ini

    T_atm_0 = CC.Fields.FieldVector(data=field_atm)
    T_oce_0 = CC.Fields.FieldVector(data=field_oce)
    T_ice_0 = CC.Fields.FieldVector(data=field_ice)
    h_ice_0 = CC.Fields.FieldVector(data=field_h_I)
    stable_range = initial_value_range([T_atm_0, T_oce_0, T_ice_0])

    if p.boundary_mapping == "cit"
        time_points = CC.Domains.IntervalDomain(
            CC.Geometry.ZPoint(stepping.timerange[1]),
            CC.Geometry.ZPoint(stepping.timerange[2]);
            boundary_names=(:start, :end),
        )
        time_mesh = CC.Meshes.IntervalMesh(
            time_points,
            nelems=Int(stepping.timerange[2] / stepping.Δt_min),
        )
        time_space = CC.Spaces.FaceFiniteDifferenceSpace(device, time_mesh)
        T_O = T_oce_0[end] .* CC.Fields.ones(time_space)
        T_A = T_atm_0[1] .* CC.Fields.ones(time_space)
        T_Is = T_ice_0[1] .* CC.Fields.ones(time_space)
    else
        T_O = T_oce_0[end] .* CC.Fields.ones(boundary_space)
        T_A = T_atm_0[1] .* CC.Fields.ones(boundary_space)
        T_Is = T_ice_0[1] .* CC.Fields.ones(point_space_ice)
    end

    parameter_dict = Dict(key => getfield(p, key) for key ∈ fieldnames(SimulationParameters))
    atmos_cache = (; parameter_dict..., T_O=T_O, T_Is=T_Is, stable_range=stable_range)
    ocean_cache = (; parameter_dict..., T_A=T_A, stable_range=stable_range)
    atmos_sim = atmos_init(stepping, T_atm_0, center_space_atm, atmos_cache)
    ocean_sim = ocean_init(stepping, T_oce_0, center_space_oce, ocean_cache)
    ice_cache = (; parameter_dict..., T_A=T_A, T_O=T_O, stable_range=stable_range)
    ice_sim = ice_init(stepping, h_ice_0, point_space_ice, ice_cache)

    comms_ctx = Utilities.get_comms_context(Dict("device" => "auto"))
    output_dir = "output"
    mkpath(output_dir)
    dir_paths = (
        output=output_dir,
        artifacts=output_dir,
        regrid=output_dir,
        checkpoints=output_dir,
    )

    start_date = "19790301"
    date = Dates.DateTime(start_date, Dates.dateformat"yyyymmdd")
    dates = (;
        date=[date],
        date0=[date],
        date1=[Dates.firstdayofmonth(date)],
        new_month=[false],
    )


    coupler_field_names = (:T_atm_sfc, :T_oce_sfc, :T_ice)
    coupler_fields = NamedTuple{coupler_field_names}(
        ntuple(i -> CC.Fields.zeros(boundary_space), length(coupler_field_names)),
    )

    model_sims = (atmos_sim=atmos_sim, ocean_sim=ocean_sim, ice_sim=ice_sim)


    cs = Interfacer.CoupledSimulation{Float64}(
        comms_ctx,
        dates,
        boundary_space,
        coupler_fields,
        nothing, # conservation checks
        stepping.timerange,
        stepping.Δt_coupler,
        model_sims,
        (;), # callbacks
        dir_paths,
        FluxCalculator.PartitionedStateFluxes(),
        nothing, # thermo_params
        nothing, # amip_diags_handler
    )
    return cs
end
