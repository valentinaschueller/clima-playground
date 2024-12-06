include("components/atmosphere.jl")
include("components/ocean.jl")
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


function solve_coupler!(cs)
    (; Δt_cpl, tspan, comms_ctx) = cs

    cs.dates.date[1] = TimeManager.current_date(cs, tspan[begin])
    TimeManager.trigger_callback!(cs, cs.callbacks.checkpoint)

    @info("Starting coupling loop")

    for t in ((tspan[begin]+Δt_cpl):Δt_cpl:tspan[end])
        @info(cs.dates.date[1])

        ClimaComms.barrier(comms_ctx)

        FieldExchanger.step_model_sims!(cs.model_sims, t)

        update_model_sims!(cs.model_sims.atmos_sim, cs.model_sims.ocean_sim)

        cs.dates.date[1] = TimeManager.current_date(cs, t)
        TimeManager.trigger_callback!(cs, cs.callbacks.checkpoint)
    end
end

function update_model_sims!(atmos_sim::HeatEquationAtmos, ocean_sim::HeatEquationOcean)
    # Get ocean temp and use it to update atmosphere model
    ocean_T = get_field(ocean_sim, Val(:T_oce_sfc))
    update_field!(atmos_sim, Val(:T_oce_sfc), Float64(ocean_T))

    # Get atmos temp and use it to update ocean model
    atmos_T = get_field(atmos_sim, Val(:T_atm_sfc))
    update_field!(ocean_sim, Val(:T_atm_sfc), atmos_T)
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
        T_atm_ini=Float64(265),   # initial temperature [K]
        T_oce_ini=Float64(271),   # initial temperature [K]
        a_i=Float64(0),           # ice area fraction [0-1]
    )

    # initialize models
    domain_atm = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(0),
        CC.Geometry.ZPoint{Float64}(parameters.h_atm);
        boundary_names=(:bottom, :top),
    )
    context = CC.ClimaComms.context()
    mesh_atm = CC.Meshes.IntervalMesh(domain_atm, nelems=parameters.n_atm)
    device = CC.ClimaComms.device(context)
    center_space_atm = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_atm)

    domain_oce = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(-parameters.h_oce),
        CC.Geometry.ZPoint{Float64}(0);
        boundary_names=(:bottom, :top),
    )
    context = CC.ClimaComms.context()
    mesh_oce = CC.Meshes.IntervalMesh(domain_oce, nelems=parameters.n_oce)
    device = CC.ClimaComms.device(context)
    center_space_oce = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_oce)

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
    )

    T_atm_0 = CC.Fields.FieldVector(
        atm=CC.Fields.ones(Float64, center_space_atm) .* parameters.T_atm_ini,
    )
    T_oce_0 = CC.Fields.FieldVector(
        oce=CC.Fields.ones(Float64, center_space_oce) .* parameters.T_oce_ini,
    )

    atmos_cache = (; parameters..., T_sfc=parameters.T_oce_ini .* CC.Fields.ones(boundary_space))
    atmos_sim = atmos_init(stepping, T_atm_0, center_space_atm, atmos_cache)
    ocean_cache = (; parameters..., T_air=parameters.T_atm_ini .* CC.Fields.ones(boundary_space), area_fraction=1 - parameters.a_i)
    ocean_sim = ocean_init(stepping, T_oce_0, center_space_oce, ocean_cache)

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

    checkpoint_cb = TimeManager.HourlyCallback(
        dt=Float64(1),
        func=Checkpointer.checkpoint_sims,
        ref_date=[dates.date[1]],
        active=true,
    )
    callbacks = (; checkpoint=checkpoint_cb)

    coupler_field_names = (
        :T_atm_sfc,
        :T_oce_sfc,
    )
    coupler_fields = NamedTuple{coupler_field_names}(
        ntuple(i -> CC.Fields.zeros(boundary_space), length(coupler_field_names)),
    )
    model_sims = (atmos_sim=atmos_sim, ocean_sim=ocean_sim)


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
        callbacks,
        dir_paths,
        FluxCalculator.PartitionedStateFluxes(),
        nothing, # thermo_params
        nothing, # amip_diags_handler
    )

    solve_coupler!(cs)

end;
