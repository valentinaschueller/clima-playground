import Dates
import ClimaAtmos as CA
import ClimaComms
import ClimaCoupler: Interfacer, ConservationChecker, FieldExchanger, TimeManager, FluxCalculator
import YAML

export coupled_run, atmos_only

include(joinpath(pkgdir(CA), "post_processing", "ci_plots.jl"))
include(
        joinpath(
            pkgdir(CA),
            "reproducibility_tests",
            "reproducibility_utils.jl",
        ),
    )

function get_simple_coupled_sim(config::CA.AtmosConfig, config_file::String)
    output_dir = "output"
    rm(output_dir, recursive=true, force=true)
    mkpath(output_dir)
    dir_paths = (
        output_dir_root=output_dir,
        artifacts=output_dir,
        regrid=output_dir,
        checkpoints_dir=output_dir,
    )

    atmos_sim = ClimaAtmosSimulation(config)
    thermo_params = get_thermo_params(atmos_sim)

    config = YAML.load_file(config_file)
    FT = config["FLOAT_TYPE"] == "Float64" ? Float64 : Float32

    sfc_sim = simple_surface_init(config, output_dir)

    boundary_space = get_surface_space(atmos_sim)

    start_date = Dates.DateTime("19790301", Dates.dateformat"yyyymmdd")

    model_sims = (atmos_sim=atmos_sim, ice_sim=sfc_sim)

    coupler_field_names = Interfacer.default_coupler_fields()
    for sim in model_sims
        Interfacer.add_coupler_fields!(coupler_field_names, sim)
    end
    coupler_fields = Interfacer.init_coupler_fields(FT, coupler_field_names, boundary_space)

    tspan = (0.0, CA.time_to_seconds(config["t_end"]))
    cs = Interfacer.CoupledSimulation{FT}(
        Ref(start_date),
        coupler_fields,
        nothing, # conservation checks
        tspan,
        CA.time_to_seconds(config["dt"]),
        Ref(tspan[begin]),
        Ref(-1),
        model_sims,
        (;), # callbacks
        dir_paths,
        thermo_params, # thermo_params
        nothing, # diagnostic_handlers
    )
    return cs
end

function step!(cs::Interfacer.CoupledSimulation)
    # Update the current time
    cs.t[] += cs.Δt_cpl

    ## compute global energy and water conservation checks
    ## (only for slabplanet if tracking conservation is enabled)
    ConservationChecker.check_conservation!(cs)

    ## step component model simulations sequentially for one coupling timestep (Δt_cpl)
    FieldExchanger.step_model_sims!(cs)

    ## update the surface fractions for surface models
    FieldExchanger.update_surface_fractions!(cs)

    ## exchange all non-turbulent flux fields between models
    FieldExchanger.exchange!(cs)

    ## calculate turbulent fluxes in the coupler and update the model simulations with them
    FluxCalculator.turbulent_fluxes!(cs)

    ## Maybe call the callbacks
    TimeManager.callbacks!(cs)

    # Compute coupler diagnostics
    # CD.orchestrate_diagnostics(cs)
    return nothing
end

function run!(
    cs::Interfacer.CoupledSimulation;
    precompile = (cs.tspan[end] > 2 * cs.Δt_cpl + cs.tspan[begin]),
)
    FieldExchanger.update_surface_fractions!(cs)
    
    #=
    ## Initialize Component Model Exchange

    The concrete steps for proper initialization are:
    =#

    # 1. Import static fields into the coupler fields
    FieldExchanger.import_static_fields!(cs.fields, cs.model_sims)

    # 2. Import atmospheric and surface fields into the coupler fields,
    #  then broadcast them back out to all components.
    FieldExchanger.exchange!(cs)

    # 3. Update any fields in the model caches that can only be filled after the initial exchange.
    FieldExchanger.set_caches!(cs)

    # 4. Calculate and update turbulent fluxes for each surface model,
    #  and save the weighted average in coupler fields
    FluxCalculator.turbulent_fluxes!(cs)

    # 4. Compute any ocean-sea ice fluxes
    FluxCalculator.ocean_seaice_fluxes!(cs)

    ## Precompilation of Coupling Loop
    # Here we run the entire coupled simulation for two timesteps to precompile several
    # functions for more accurate timing of the overall simulation.
    precompile && (step!(cs); step!(cs))

    ## Run garbage collection before solving for more accurate memory comparison to ClimaAtmos
    GC.gc()

    @info "Starting coupling loop"
    walltime = ClimaComms.@elapsed ClimaComms.device(cs) begin
        s = CA.@timed_str begin
            while cs.t[] <= cs.tspan[end]
                step!(cs)
            end
        end
    end
    @info "Simulation took $(walltime) seconds"

    # Close all diagnostics file writers
    isnothing(cs.diags_handler) ||
        foreach(diag -> close(diag.output_writer), cs.diags_handler.scheduled_diagnostics)
    foreach(Interfacer.close_output_writers, cs.model_sims)

    return nothing
end

function coupled_run(;
    config_file::String,
    reference_job_id::String,
)
    redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))

    config = CA.AtmosConfig(config_file)
    cs = get_simple_coupled_sim(config, config_file)

    run!(cs)

    output_dir = cs.dir_paths.output_dir_root
    @info "Plotting"
    make_plots(Val(Symbol(reference_job_id)), output_dir)
    @info "Plotting done"

    return cs
end

function atmos_only(;
    config_file::String,
    reference_job_id::String,
)
    redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))

    config = CA.AtmosConfig(config_file)
    simulation = CA.get_simulation(config)
    CA.solve_atmos!(simulation)

    output_dir = simulation.output_dir
    @info "Plotting"
    make_plots(Val(Symbol(reference_job_id)), output_dir)
    @info "Plotting done"
end
