include("components/atmosphere_with_ice.jl")
include("components/ocean_with_ice.jl")
include("components/ice.jl")
import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler:
    Checkpointer,
    FieldExchanger,
    FluxCalculator,
    Interfacer,
    TimeManager,
    Utilities
# This is not good!!! Or maybe it is, if the otherone cannot checkpoint at the individual steps, which i would want.


function reset_time!(cs::Interfacer.CoupledSimulation, t)
    for sim in cs.model_sims
        sim.integrator.t = t
    end
end

# function restart_sims!(cs::Interfacer.CoupledSimulation)
#     @info "Reading checkpoint!"
#     for sim in cs.model_sims
#         if Checkpointer.get_model_prog_state(sim) !== nothing
#             time = Dates.datetime2epochms(cs.dates.date[1])
#             t0 = Dates.datetime2epochms(cs.dates.date0[1])
#             Checkpointer.restart_model_state!(sim, cs.comms_ctx, Int((time - t0) / 1e3), input_dir=cs.dirs.artifacts)
#         end
#     end
# end

function solve_coupler!(cs::Interfacer.CoupledSimulation, max_iters)
    (; Δt_cpl, tspan) = cs

    t_range = ((tspan[begin]+Δt_cpl):Δt_cpl:tspan[end])
    bound_temps = [fill(get_field(cs.model_sims.atmos_sim, Val(:T_atm_sfc)), length(t_range) + 1), fill(get_field(cs.model_sims.ocean_sim, Val(:T_oce_sfc)), length(t_range) + 1), fill(get_field(cs.model_sims.ice_sim, Val(:T_ice)), length(t_range) + 1)]

    @info("Starting coupling loop")
    iter = 1
    time_vec = [0]

    while iter <= max_iters
        @info("Current iter: $(iter)")
        cs.dates.date[1] = TimeManager.current_date(cs, tspan[begin])

        time = Dates.datetime2epochms(cs.dates.date[1])
        t0 = Dates.datetime2epochms(cs.dates.date0[1])
        Checkpointer.checkpoint_model_state(cs.model_sims.atmos_sim, cs.comms_ctx, Int((time - t0) / 1e3), output_dir=cs.dirs.artifacts)
        Checkpointer.checkpoint_model_state(cs.model_sims.ocean_sim, cs.comms_ctx, Int((time - t0) / 1e3), output_dir=cs.dirs.artifacts)

        for (i, t) in enumerate(t_range)

            update_field!(cs.model_sims.atmos_sim, Val(:T_oce_sfc), bound_temps[2][i], Val(:T_ice), bound_temps[3][i])

            Interfacer.step!(cs.model_sims.atmos_sim, t)

            atmos_T = get_field(cs.model_sims.atmos_sim, Val(:T_atm_sfc))

            bound_temps[1][i+1] = atmos_T

            cs.dates.date[1] = TimeManager.current_date(cs, t)

            time = Dates.datetime2epochms(cs.dates.date[1])
            t0 = Dates.datetime2epochms(cs.dates.date0[1])
            push!(time_vec, Int((time - t0) / 1e3))
            Checkpointer.checkpoint_model_state(cs.model_sims.atmos_sim, cs.comms_ctx, Int((time - t0) / 1e3), output_dir=cs.dirs.artifacts)
        end

        for (i, t) in enumerate(t_range)
            update_field!(cs.model_sims.ocean_sim, Val(:T_atm_sfc), bound_temps[1][i], Val(:T_ice), bound_temps[3][i])

            Interfacer.step!(cs.model_sims.ocean_sim, t)

            ocean_T = get_field(cs.model_sims.ocean_sim, Val(:T_oce_sfc))

            println(ocean_T)

            bound_temps[2][i+1] = ocean_T

            cs.dates.date[1] = TimeManager.current_date(cs, t)

            Checkpointer.checkpoint_model_state(cs.model_sims.ocean_sim, cs.comms_ctx, time_vec[i+1], output_dir=cs.dirs.artifacts)
        end

        for i in 1:1:length(t_range)+1
            t = time_vec[i]
            original_file = joinpath(cs.dirs.artifacts, "checkpoint", "checkpoint_" * Interfacer.name(cs.model_sims.ocean_sim) * "_$t.hdf5")
            new_file = joinpath(cs.dirs.artifacts, "checkpoint", "checkpoint_" * Interfacer.name(cs.model_sims.ocean_sim) * "_$iter" * "_$t.hdf5")
            print(i)
            mv(original_file, new_file, force=true)
            println("hello")
            original_file = joinpath(cs.dirs.artifacts, "checkpoint", "checkpoint_" * Interfacer.name(cs.model_sims.atmos_sim) * "_$t.hdf5")
            new_file = joinpath(cs.dirs.artifacts, "checkpoint", "checkpoint_" * Interfacer.name(cs.model_sims.atmos_sim) * "_$iter" * "_$t.hdf5")
            mv(original_file, new_file, force=true)
        end

        iter += 1
        if iter <= max_iters
            # restart_sims!(cs)
            reset_time!(cs, tspan[begin])
        end
    end
end

function coupled_heat_equations()
    parameters = (
        h_atm=Float64(200),   # depth [m]
        h_oce=Float64(50),    # depth [m]
        n_atm=200,
        n_oce=50,
        k_atm=Float64(0.02364),
        k_oce=Float64(0.57),
        c_atm=Float64(1e-3),  # specific heat [J / kg / K]
        c_oce=Float64(800),   # specific heat [J / kg / K]
        ρ_atm=Float64(1),     # density [kg / m3]
        ρ_oce=Float64(1000),  # density [kg / m3]
        C_AO=Float64(1e-5),
        C_AI=Float64(1e-5),
        C_OI=Float64(1e-5),
        T_atm_ini=Float64(268),   # initial temperature [K]
        T_oce_ini=Float64(264),   # initial temperature [K]
        T_ice_ini=Float64(260),       # temperature [K]
        a_i=Float64(0),           # ice area fraction [0-1]
    )

    context = CC.ClimaComms.context()
    device = CC.ClimaComms.device(context)

    # initialize models
    domain_atm = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(0),
        CC.Geometry.ZPoint{Float64}(parameters.h_atm);
        boundary_names=(:bottom, :top),
    )
    mesh_atm = CC.Meshes.IntervalMesh(domain_atm, nelems=parameters.n_atm)
    center_space_atm = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_atm)

    domain_oce = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(-parameters.h_oce),
        CC.Geometry.ZPoint{Float64}(0);
        boundary_names=(:bottom, :top),
    )
    mesh_oce = CC.Meshes.IntervalMesh(domain_oce, nelems=parameters.n_oce)
    center_space_oce = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_oce)

    # domain_ice = CC.Domains.IntervalDomain(
    #     CC.Geometry.ZPoint{Float64}(-0.1),
    #     CC.Geometry.ZPoint{Float64}(0.1);
    #     boundary_names=(:bottom, :top),
    # )
    # mesh_ice = CC.Meshes.IntervalMesh(domain_ice, nelems=1)
    coord = CC.Geometry.ZPoint{Float64}(0.0)
    point_space_ice = CC.Spaces.PointSpace(context, coord)

    atmos_facespace = CC.Spaces.FaceFiniteDifferenceSpace(center_space_atm)
    boundary_space = CC.Spaces.level(
        atmos_facespace,
        CC.Utilities.PlusHalf(CC.Spaces.nlevels(atmos_facespace) - 1),
    )

    stepping = (;
        Δt_min=Float64(1.0),
        timerange=(Float64(0.0), Float64(3600.0)),
        Δt_coupler=Float64(100.0),
        odesolver=CTS.ExplicitAlgorithm(CTS.RK4()),
        nsteps_atm=50,
        nsteps_oce=1,
        nsteps_ice=1,
    )

    T_atm_0 = CC.Fields.FieldVector(
        atm=CC.Fields.ones(Float64, center_space_atm) .* parameters.T_atm_ini,
    )
    T_oce_0 = CC.Fields.FieldVector(
        oce=CC.Fields.ones(Float64, center_space_oce) .* parameters.T_oce_ini,
    )
    T_ice_0 = CC.Fields.FieldVector(
        ice=CC.Fields.ones(Float64, point_space_ice) .* parameters.T_ice_ini,
    )

    atmos_cache = (; parameters..., T_sfc=parameters.T_oce_ini .* CC.Fields.ones(boundary_space), T_ice=T_ice_0)
    atmos_sim = atmos_init(stepping, T_atm_0, center_space_atm, atmos_cache)
    ocean_cache = (; parameters..., T_air=parameters.T_atm_ini .* CC.Fields.ones(boundary_space), T_ice=T_ice_0)
    ocean_sim = ocean_init(stepping, T_oce_0, center_space_oce, ocean_cache)
    ice_cache = (; parameters...)
    ice_sim = ice_init(stepping, T_ice_0, point_space_ice, ice_cache)

    comms_ctx = Utilities.get_comms_context(Dict("device" => "auto"))
    dir_paths = (output=".", artifacts=".", regrid=".")

    start_date = "19790301"
    date = Dates.DateTime(start_date, Dates.dateformat"yyyymmdd")
    dates = (;
        date=[date],
        date0=[date],
        date1=[Dates.firstdayofmonth(date)],
        new_month=[false],
    )


    coupler_field_names = (
        :T_atm_sfc,
        :T_oce_sfc,
        :T_ice,
    )
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
        (;), # mode_specifics
        (;),
        dir_paths,
        FluxCalculator.PartitionedStateFluxes(),
        nothing, # thermo_params
        nothing, # amip_diags_handler
    )

    solve_coupler!(cs, 3)

end;
