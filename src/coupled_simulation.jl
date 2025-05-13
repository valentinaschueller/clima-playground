import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities


function get_coupled_sim(p::SimulationParameters)
    context = CC.ClimaComms.context()
    device = CC.ClimaComms.device(context)

    domain_A = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(p.z_A0),
        CC.Geometry.ZPoint{Float64}(p.h_A);
        boundary_names=(:bottom, :top),
    )
    mesh_A = CC.Meshes.IntervalMesh(domain_A, nelems=p.n_A)
    center_space_atm = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_A)

    domain_O = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(-p.h_O),
        CC.Geometry.ZPoint{Float64}(-p.z_O0);
        boundary_names=(:bottom, :top),
    )
    mesh_O = CC.Meshes.IntervalMesh(domain_O, nelems=p.n_O)
    center_space_oce = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_O)

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
        nsteps_ice=1,
    )

    field_atm = CC.Fields.ones(Float64, center_space_atm) .* p.T_A_ini
    field_oce = CC.Fields.ones(Float64, center_space_oce) .* p.T_O_ini
    field_ice = CC.Fields.ones(Float64, point_space_ice) .* p.T_I_ini

    T_atm_0 = CC.Fields.FieldVector(data=field_atm)
    T_oce_0 = CC.Fields.FieldVector(data=field_oce)
    T_ice_0 = CC.Fields.FieldVector(data=field_ice)

    if p.boundary_mapping == "cit"
        time_points = CC.Domains.IntervalDomain(
            CC.Geometry.ZPoint{Float64}(stepping.timerange[1] - (stepping.Δt_min / 2)),
            CC.Geometry.ZPoint{Float64}(stepping.timerange[2] - (stepping.Δt_min / 2));
            boundary_names=(:start, :end),
        )
        mesh_time = CC.Meshes.IntervalMesh(
            time_points,
            nelems=Int(stepping.timerange[2] / stepping.Δt_min + 1),
        )
        space_time = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_time)
        T_sfc = p.T_O_ini .* CC.Fields.ones(space_time)
        T_air = p.T_A_ini .* CC.Fields.ones(space_time)
    else
        T_sfc = p.T_O_ini .* CC.Fields.ones(boundary_space)
        T_air = p.T_A_ini .* CC.Fields.ones(boundary_space)
    end

    parameter_dict = Dict(key => getfield(p, key) for key ∈ fieldnames(SimulationParameters))
    atmos_cache = (; parameter_dict..., T_sfc=T_sfc, T_ice=T_ice_0)
    ocean_cache = (; parameter_dict..., T_air=T_air, T_ice=T_ice_0)
    atmos_sim = atmos_init(stepping, T_atm_0, center_space_atm, atmos_cache)
    ocean_sim = ocean_init(stepping, T_oce_0, center_space_oce, ocean_cache)
    ice_cache = (; parameter_dict...)
    ice_sim = ice_init(stepping, T_ice_0, point_space_ice, ice_cache)

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
