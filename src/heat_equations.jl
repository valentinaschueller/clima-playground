include("components/atmosphere.jl")
include("components/ocean.jl")
import Dates
import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler
import ClimaCoupler:
    ConservationChecker,
    Checkpointer,
    Diagnostics,
    FieldExchanger,
    FluxCalculator,
    Interfacer,
    Regridder,
    TimeManager,
    Utilities


function initialize_component_model_exchange!(cs::Interfacer.CoupledSimulation)
    # 1.surface density (`ρ_sfc`): calculated by the coupler by adiabatically extrapolating atmospheric thermal state to the surface.
    ## For this we need to import surface and atmospheric fields. The model sims are then updated with the new surface density.
    FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes)
    FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes)
    FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

    # 2.surface vapor specific humidity (`q_sfc`): step surface models with the new surface density to calculate their respective `q_sfc` internally
    Interfacer.step!(ocean_sim, Δt_cpl)

    # 3.turbulent fluxes
    ## import the new surface properties into the coupler (note the atmos state was also imported in step 3.)
    FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes) # i.e. T_sfc, albedo, z0, beta, q_sfc
    ## calculate turbulent fluxes inside the atmos cache based on the combined surface state in each grid box
    FluxCalculator.combined_turbulent_fluxes!(cs.model_sims, cs.fields, cs.turbulent_fluxes) # this updates the atmos thermo state, sfc_ts

    # 4.reinitialize models + radiative flux: prognostic states and time are set to their initial conditions.
    FieldExchanger.reinit_model_sims!(cs.model_sims)

    # 5.update all fluxes: coupler re-imports updated atmos fluxes
    FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes)
    FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)
end

function solve_coupler!(cs)
    (; Δt_cpl, tspan, comms_ctx) = cs

    @info("Starting coupling loop")
    ## step in time
    for t in ((tspan[begin] + Δt_cpl):Δt_cpl:tspan[end])

        cs.dates.date[1] = TimeManager.current_date(cs, t)

        ## print date on the first of month
        cs.dates.date[1] >= cs.dates.date1[1] && @info(cs.dates.date[1])

        ClimaComms.barrier(comms_ctx)

        ## run component models sequentially for one coupling timestep (Δt_cpl)
        FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

        ## step sims
        FieldExchanger.step_model_sims!(cs.model_sims, t)

        ## exchange combined fields and (if specified) calculate fluxes using combined states
        FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes) # i.e. T_sfc, surface_albedo, z0, beta
        FluxCalculator.combined_turbulent_fluxes!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

        FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes) # radiative and/or turbulent

        ## callback to update the fist day of month if needed
        TimeManager.trigger_callback!(cs, cs.callbacks.update_firstdayofmonth!)

        ## callback to checkpoint model state
        TimeManager.trigger_callback!(cs, cs.callbacks.checkpoint)

    end

    return nothing
end


function coupled_heat_equations()
    parameters = (
        h_atm = Float64(1.0),   # depth [m]
        h_oce = Float64(20),    # depth [m]
        n_atm = 15,
        n_oce = 15,
        k_atm = Float64(1e-3),
        k_oce = Float64(1e-3),
        c_atm = Float64(1e-3),  # specific heat [J / kg / K]
        c_oce = Float64(800),   # specific heat [J / kg / K]
        ρ_atm = Float64(1),     # density [kg / m3]
        ρ_oce = Float64(1500),  # density [kg / m3]
        C_AO = Float64(1e-5),
        C_AI = Float64(1e-5),
        C_OI = Float64(1e-5),
        T_atm_ini = Float64(265),   # initial temperature [K]
        T_oce_ini = Float64(271),   # initial temperature [K]
        a_i = Float64(0),           # ice area fraction [0-1]
    )

    # initialize models
    domain_atm = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(0),
        CC.Geometry.ZPoint{Float64}(parameters.h_atm);
        boundary_names = (:bottom, :top),
    );
    context = CC.ClimaComms.context()
    mesh_atm = CC.Meshes.IntervalMesh(domain_atm, nelems = parameters.n_atm); # struct, allocates face boundaries to 5,6: atmos
    device = CC.ClimaComms.device(context)
    center_space_atm = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_atm) 

    domain_oce = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(-parameters.h_oce),
        CC.Geometry.ZPoint{Float64}(0);
        boundary_names = (:bottom, :top),
    );
    context = CC.ClimaComms.context()
    mesh_oce = CC.Meshes.IntervalMesh(domain_oce, nelems = parameters.n_oce); # struct, allocates face boundaries to 5,6: atmos
    device = CC.ClimaComms.device(context)
    center_space_oce = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_oce) 

    stepping = (;
        Δt_min = Float64(0.02),
        timerange = (Float64(0.0), Float64(6.0)),
        Δt_coupler = Float64(0.1),
        odesolver = CTS.ExplicitAlgorithm(CTS.RK4()),
        nsteps_atm = 8, # number of timesteps of atm per coupling cycle
        nsteps_oce = 1, # number of timesteps of lnd per coupling cycle
    );

    T_atm_0 = CC.Fields.FieldVector(atm = CC.Fields.ones(Float64, center_space_atm) .* parameters.T_atm_ini);
    T_oce_0 = CC.Fields.FieldVector(oce = CC.Fields.ones(Float64, center_space_oce) .* parameters.T_oce_ini);

    ## atmos copies of coupler variables
    atmos_sim = atmos_init(stepping, T_atm_0, center_space_atm, parameters)
    ocean_sim = ocean_init(stepping, T_oce_0, center_space_oce, parameters)

    comms_ctx = Utilities.get_comms_context(Dict("device" => "auto"))

    start_date = "19790301"
    date0 = date = Dates.DateTime(start_date, Dates.dateformat"yyyymmdd")
    dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)], new_month = [false])

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
    coupler_fields =
        NamedTuple{coupler_field_names}(ntuple(i -> CC.Fields.zeros(boundary_space), length(coupler_field_names)))

    model_sims = (atmos_sim = atmos_sim, ocean_sim = ocean_sim)
    

    cs = Interfacer.CoupledSimulation{Float64}(
        comms_ctx,
        dates,
        CC.Spaces.horizontal_space(atmos_sim.domain.face_space),
        coupler_fields,
        nothing, # conservation checks
        stepping.timerange,
        stepping.Δt_coupler,
        model_sims,
        (;), # mode_specifics
        callbacks,
        dir_paths,
        turbulent_fluxes,
        thermo_params,
        nothing, # amip_diags_handler
    );

    initialize_component_model_exchange!(cs)

    integ_atm, integ_oce = solve_coupler!(cs)

    # postprocessing
    sol_atm, sol_lnd = integ_atm.sol, integ_oce.sol;
end;


