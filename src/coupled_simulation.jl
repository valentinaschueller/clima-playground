import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities

"""
Creates a CoupledSimulation.

**Arguments:**

-`physical_values::Dict`: Can be defined using `define_realistic_vals()`. 

"""
function get_coupled_sim(parameters)
    context = CC.ClimaComms.context()
    device = CC.ClimaComms.device(context)

    domain_atm = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(parameters[:z_0numA]),
        CC.Geometry.ZPoint{Float64}(parameters[:h_atm]);
        boundary_names=(:bottom, :top),
    )
    mesh_atm = CC.Meshes.IntervalMesh(domain_atm, nelems=parameters[:n_atm])
    center_space_atm = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_atm)

    domain_oce = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(-parameters[:h_oce]),
        CC.Geometry.ZPoint{Float64}(-parameters[:z_0numO]);
        boundary_names=(:bottom, :top),
    )
    mesh_oce = CC.Meshes.IntervalMesh(domain_oce, nelems=parameters[:n_oce])
    center_space_oce = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_oce)

    coord = CC.Geometry.ZPoint{Float64}(0.0)
    point_space_ice = CC.Spaces.PointSpace(context, coord)

    atmos_facespace = CC.Spaces.FaceFiniteDifferenceSpace(center_space_atm)
    boundary_space = CC.Spaces.level(
        atmos_facespace,
        CC.Utilities.PlusHalf(CC.Spaces.nlevels(atmos_facespace) - 1),
    )

    stepping = (;
        Δt_min=Float64(parameters[:Δt_min]),
        timerange=(Float64(0.0), Float64(parameters[:t_max])),
        Δt_coupler=Float64(parameters[:Δt_cpl]),
        odesolver=CTS.ExplicitAlgorithm(CTS.RK4()),
        nsteps_atm=parameters[:n_t_atm],
        nsteps_oce=parameters[:n_t_oce],
        nsteps_ice=1,
    )
    if parameters[:sin_field_atm]
        coord_field_atm = map(x -> x.z, CC.Fields.coordinate_field(center_space_atm))
        field_atm =
            parameters[:T_atm_ini] .* (
                1 .-
                sin.(
                    (coord_field_atm .- parent(coord_field_atm)[1]) ./
                    (parent(coord_field_atm)[end] .- parent(coord_field_atm)[1]) .*
                    (π / 50)
                )
            )
    else
        field_atm = CC.Fields.ones(Float64, center_space_atm) .* parameters[:T_atm_ini]
    end

    if parameters[:sin_field_oce]
        coord_field_oce = map(x -> x.z, CC.Fields.coordinate_field(center_space_oce))
        field_oce =
            1 .+
            sin.(
                (coord_field_oce .- parent(coord_field_oce)[1]) ./
                (parent(coord_field_oce)[end] .- parent(coord_field_oce)[1]) .* (π / 50)
            )
        field_oce = parameters[:T_oce_ini] .* (field_oce .- (parent(field_oce)[end] - 1))

    else
        field_oce = CC.Fields.ones(Float64, center_space_oce) .* parameters[:T_oce_ini]
    end

    T_atm_0 = CC.Fields.FieldVector(atm=field_atm)
    T_oce_0 = CC.Fields.FieldVector(oce=field_oce)
    T_ice_0 = CC.Fields.FieldVector(
        ice=CC.Fields.ones(Float64, point_space_ice) .* parameters[:T_ice_ini],
    )

    if parameters[:boundary_mapping] == "cit"
        time_points_oce_domain = CC.Domains.IntervalDomain(
            CC.Geometry.ZPoint{Float64}(stepping.timerange[1] - (stepping.Δt_min / 2)),
            CC.Geometry.ZPoint{Float64}(stepping.timerange[2] - (stepping.Δt_min / 2));
            boundary_names=(:start, :end),
        )
        time_points_atm_domain = CC.Domains.IntervalDomain(
            CC.Geometry.ZPoint{Float64}(stepping.timerange[1] - (stepping.Δt_min / 2)),
            CC.Geometry.ZPoint{Float64}(stepping.timerange[2] - (stepping.Δt_min / 2));
            boundary_names=(:start, :end),
        )

        mesh_time_oce = CC.Meshes.IntervalMesh(
            time_points_oce_domain,
            nelems=Int(stepping.timerange[2] / stepping.Δt_min + 1),
        )
        mesh_time_atm = CC.Meshes.IntervalMesh(
            time_points_atm_domain,
            nelems=Int(stepping.timerange[2] / stepping.Δt_min + 1),
        )
        space_time_oce = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_time_oce)
        space_time_atm = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_time_atm)
        T_sfc = parameters[:T_oce_ini] .* CC.Fields.ones(space_time_oce)
        T_air = parameters[:T_atm_ini] .* CC.Fields.ones(space_time_atm)
    else
        T_sfc = parameters[:T_oce_ini] .* CC.Fields.ones(boundary_space)
        T_air = parameters[:T_atm_ini] .* CC.Fields.ones(boundary_space)
    end

    atmos_cache = (; parameters..., T_sfc=T_sfc, T_ice=T_ice_0)
    ocean_cache = (; parameters..., T_air=T_air, T_ice=T_ice_0)
    atmos_sim = atmos_init(stepping, T_atm_0, center_space_atm, atmos_cache)
    ocean_sim = ocean_init(stepping, T_oce_0, center_space_oce, ocean_cache)
    ice_cache = (; parameters...)
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
