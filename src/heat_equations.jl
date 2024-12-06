include("components/atmosphere.jl")
include("components/ocean.jl")
import Dates
import SciMLBase
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler
import ClimaCoupler:
    ConservationChecker,
    Checkpointer,
    FieldExchanger,
    FluxCalculator,
    Interfacer,
    Regridder,
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

        ## run component models sequentially for one coupling timestep (Δt_cpl)
        FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

        ## step sims
        FieldExchanger.step_model_sims!(cs.model_sims, t)

        ## exchange combined fields and (if specified) calculate fluxes using combined states
        FieldExchanger.import_combined_surface_fields!(
            cs.fields,
            cs.model_sims,
            cs.turbulent_fluxes,
        ) # i.e. T_sfc, surface_albedo, z0, beta
        FluxCalculator.combined_turbulent_fluxes!(
            cs.model_sims,
            cs.fields,
            cs.turbulent_fluxes,
        )

        FieldExchanger.import_atmos_fields!(
            cs.fields,
            cs.model_sims,
            cs.boundary_space,
            cs.turbulent_fluxes,
        ) # radiative and/or turbulent

        cs.dates.date[1] = TimeManager.current_date(cs, t)
        TimeManager.trigger_callback!(cs, cs.callbacks.checkpoint)
    end
end


function coupled_heat_equations()
    parameters = (
        h_atm=Float64(1),   # depth [m]
        h_oce=Float64(20),    # depth [m]
        n_atm=15,
        n_oce=15,
        k_atm=Float64(1e-3),
        k_oce=Float64(1e-3),
        c_atm=Float64(1e-3),  # specific heat [J / kg / K]
        c_oce=Float64(800),   # specific heat [J / kg / K]
        ρ_atm=Float64(1),     # density [kg / m3]
        ρ_oce=Float64(1500),  # density [kg / m3]
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
    mesh_atm = CC.Meshes.IntervalMesh(domain_atm, nelems=parameters.n_atm) # struct, allocates face boundaries to 5,6: atmos
    device = CC.ClimaComms.device(context)
    center_space_atm = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_atm)

    domain_oce = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(-parameters.h_oce),
        CC.Geometry.ZPoint{Float64}(0);
        boundary_names=(:bottom, :top),
    )
    context = CC.ClimaComms.context()
    mesh_oce = CC.Meshes.IntervalMesh(domain_oce, nelems=parameters.n_oce) # struct, allocates face boundaries to 5,6: atmos
    device = CC.ClimaComms.device(context)
    center_space_oce = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_oce)

    stepping = (;
        Δt_min=Float64(1.0),
        timerange=(Float64(0.0), Float64(3600.0)),
        Δt_coupler=Float64(3600.0),
        odesolver=CTS.ExplicitAlgorithm(CTS.RK4()),
        nsteps_atm=1, # number of timesteps of atm per coupling cycle
        nsteps_oce=1, # number of timesteps of lnd per coupling cycle
    )

    T_atm_0 = CC.Fields.FieldVector(
        atm=CC.Fields.ones(Float64, center_space_atm) .* parameters.T_atm_ini,
    )
    T_oce_0 = CC.Fields.FieldVector(
        oce=CC.Fields.ones(Float64, center_space_oce) .* parameters.T_oce_ini,
    )

    ## atmos copies of coupler variables
    atmos_sim = atmos_init(stepping, T_atm_0, center_space_atm, parameters)
    ocean_sim = ocean_init(stepping, T_oce_0, center_space_oce, parameters)

    comms_ctx = Utilities.get_comms_context(Dict("device" => "auto"))
    dir_paths = Utilities.setup_output_dirs(
        output_dir=".",
        artifacts_dir=".",
        comms_ctx=comms_ctx,
    )

    start_date = "19790301"
    date0 = date = Dates.DateTime(start_date, Dates.dateformat"yyyymmdd")
    dates = (;
        date=[date],
        date0=[date0],
        date1=[Dates.firstdayofmonth(date0)],
        new_month=[false],
    )

    atmos_facespace = CC.Spaces.FaceFiniteDifferenceSpace(center_space_atm)
    boundary_space = CC.Spaces.level(
        atmos_facespace,
        CC.Utilities.PlusHalf(CC.Spaces.nlevels(atmos_facespace) - 1),
    )

    checkpoint_cb = TimeManager.HourlyCallback(
        dt=Float64(1),
        func=Checkpointer.checkpoint_sims,
        ref_date=[dates.date[1]],
        active=true,
    )
    callbacks = (; checkpoint=checkpoint_cb)

    coupler_field_names = (
        :T_S,
        :z0m_S,
        :z0b_S,
        :ρ_sfc,
        :q_sfc,
        :surface_direct_albedo,
        :surface_diffuse_albedo,
        :beta,
        :F_turb_energy,
        :F_turb_moisture,
        :F_turb_ρτxz,
        :F_turb_ρτyz,
        :F_radiative,
        :P_liq,
        :P_snow,
        :radiative_energy_flux_toa,
        :P_net,
        :temp1,
        :temp2,
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
