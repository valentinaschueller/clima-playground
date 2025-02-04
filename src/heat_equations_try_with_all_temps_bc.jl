include("components/atmosphere_with_ice_try_bc.jl")
include("components/ocean_with_ice_try_bc.jl")
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

function rename_file(cs::Interfacer.CoupledSimulation, iter, time, reverse=false)
    """When a file is saved, its always called the same thing.
    Had to rename it for each iteration to not overwrite"""

    original_file = joinpath(cs.dirs.artifacts, "checkpoint", "checkpoint_" * Interfacer.name(cs.model_sims.ocean_sim) * "_$time.hdf5")
    new_file = joinpath(cs.dirs.artifacts, "checkpoint", "checkpoint_" * Interfacer.name(cs.model_sims.ocean_sim) * "_$iter" * "_$time.hdf5")
    if !reverse
        mv(original_file, new_file, force=true)
    else
        mv(new_file, original_file, force=true)
    end
    original_file = joinpath(cs.dirs.artifacts, "checkpoint", "checkpoint_" * Interfacer.name(cs.model_sims.atmos_sim) * "_$time.hdf5")
    new_file = joinpath(cs.dirs.artifacts, "checkpoint", "checkpoint_" * Interfacer.name(cs.model_sims.atmos_sim) * "_$iter" * "_$time.hdf5")
    if !reverse
        mv(original_file, new_file, force=true)
    else
        mv(new_file, original_file, force=true)
    end
end

function reset_time!(cs::Interfacer.CoupledSimulation, t)
    """resets integrator time"""
    for sim in cs.model_sims
        sim.integrator.t = t
    end
end

function restart_sims!(cs::Interfacer.CoupledSimulation)
    @info "Reading checkpoint!"
    for sim in cs.model_sims
        if Checkpointer.get_model_prog_state(sim) !== nothing
            t = Dates.datetime2epochms(cs.dates.date[1])
            t0 = Dates.datetime2epochms(cs.dates.date0[1])
            time = Int((t - t0) / 1e3)

            rename_file(cs, 0, time, true)
            Checkpointer.restart_model_state!(sim, cs.comms_ctx, time, input_dir=cs.dirs.artifacts)
            rename_file(cs, 0, time)
        end
    end
    # Had to add this to restart the integrator and get access to the new temperatures in 
    # each iteration.
    for sim in cs.model_sims
        Interfacer.reinit!(sim.integrator, sim.Y_init)
    end
end

function solve_coupler!(cs::Interfacer.CoupledSimulation, tol, print_conv, plot_conv)
    (; Δt_cpl, tspan) = cs

    cs.dates.date[1] = TimeManager.current_date(cs, tspan[begin])

    @info("Starting coupling loop")

    for t in ((tspan[begin]+Δt_cpl):Δt_cpl:tspan[end])

        time = Int(t - Δt_cpl)
        @info(cs.dates.date[1])
        diff_atm = typemax(Float64)
        diff_oce = typemax(Float64)

        Checkpointer.checkpoint_sims(cs) # I had to remove nothing here
        rename_file(cs, 0, time)

        times = []
        numeric_z_range_ocean = []
        numeric_z_range_atmos = []
        times = cs.model_sims.ocean_sim.integrator.sol.t
        z_range_ocean = cs.model_sims.ocean_sim.domain.grid.topology.mesh.faces
        z_range_atmos = cs.model_sims.atmos_sim.domain.grid.topology.mesh.faces
        for i in 1:1:length(z_range_atmos)
            push!(numeric_z_range_atmos, z_range_atmos[i].z)
        end
        for i in 1:1:length(z_range_ocean)
            push!(numeric_z_range_ocean, z_range_ocean[i].z)
        end

        iter = 1
        atmos_vals_list = []
        ocean_vals_list = []
        atmos_vals = Nothing
        ocean_vals = Nothing

        while diff_oce > tol && diff_atm > tol
            @info("Current iter: $(iter)")
            # Update models
            FieldExchanger.step_model_sims!(cs.model_sims, t)

            atmos_states = copy(cs.model_sims.atmos_sim.integrator.sol.u)
            ocean_states = copy(cs.model_sims.ocean_sim.integrator.sol.u)

            # Temperature values for this iteration.
            pre_atmos_vals = atmos_vals
            pre_ocean_vals = ocean_vals
            atmos_vals = extract_matrix(atmos_states, "atm")
            ocean_vals = extract_matrix(ocean_states, "oce")
            push!(atmos_vals_list, atmos_vals)
            push!(ocean_vals_list, ocean_vals)
            if iter > 1
                diff_atm = norm(atmos_vals .- pre_atmos_vals)
                diff_oce = norm(ocean_vals .- pre_ocean_vals)
            end

            # Save temperatures
            save_temp(atmos_states, "atm", iter, numeric_z_range_atmos[1:end-1], times)
            save_temp(ocean_states, "oce", iter, numeric_z_range_ocean[2:end], times)
            # save_temp(ice_states, "ice", iter, [0], times)

            Checkpointer.checkpoint_sims(cs) # I had to remove nothing here

            rename_file(cs, iter, time)

            ice_T = get_field(cs.model_sims.ice_sim, Val(:T_ice))

            restart_sims!(cs)
            reset_time!(cs, t - Δt_cpl)
            set_coupling_fields!(cs.model_sims.atmos_sim, cs.model_sims.ocean_sim, cs.model_sims.ice_sim, atmos_states, ocean_states, ice_T)
            iter = iter + 1
        end
        cs.dates.date[1] = TimeManager.current_date(cs, t)

        if print_conv || plot_conv
            conv_fac_atm = []
            conv_fac_oce = []
            error_atm = 0
            error_oce = 0
            for i = 1:iter-1
                pre_error_atm = error_atm
                pre_error_oce = error_oce
                full_errors_atm = atmos_vals_list[i] .- atmos_vals_list[end]
                error_atm = maximum([abs(full_error_atm[1]) for full_error_atm in full_errors_atm])
                full_errors_oce = ocean_vals_list[i] .- ocean_vals_list[end]
                error_oce = maximum([abs(full_error_oce[end]) for full_error_oce in full_errors_oce])
                if i > 1
                    push!(conv_fac_atm, error_atm / pre_error_atm)
                    push!(conv_fac_oce, error_oce / pre_error_oce)
                end
            end
            if print_conv
                println("Convergence factor atmosphere: $conv_fac_atm")
                println("Convergence factor atmosphere: $conv_fac_atm")
            end
            if plot_conv
                k = 2:iter-1
                scatter(k, conv_fac_atm, label="atm", color=:blue, markersize=5, xlabel="Iteration for last used temperature", ylabel="Convergence factor")
                scatter!(k, conv_fac_oce, label="oce", color=:green, markersize=5)
                display(current())
            end # Allow for computation, plot and print of convergence factor in this script.
        end # Use mean of convergence factors over iteration? plot as a function of deltat...
    end
    # Thought: the Schwarz iteration I do here is not the same as in my analysis right? It's the "simultaneous version"

    # Todo ideas: I guess it would be interesting to run for a very long time, as we then might have heat
    # propagating further into the domain. Now there is no difference right, it has not propagated yet
    # Change with iterations. I saw something that might imply that we would increase and then decrease with increased iterations.
    # Also change in initial conditions to maybe see a larger difference in other areas of the domain.
    # Check the boundary condition at z=0 with valentina. Write down stuff. try adding an update formula for the ice.
end

function set_coupling_fields!(atmos_sim::HeatEquationAtmos, ocean_sim::HeatEquationOcean, ice_sim::ConstantIce, atmos_T, ocean_T, ice_T)
    update_field!(ice_sim, Val(:T_ice), ice_T)
    update_field!(atmos_sim, ocean_T, ice_T)
    update_field!(ocean_sim, atmos_T, ice_T)
end

function extract_matrix(field_vecs, type)
    matrix = []
    for field_vec in field_vecs
        if type == "atm"
            field = field_vec.atm
            values = parent(field)
            push!(matrix, values)
        else
            field = field_vec.oce
            values = parent(field)
            push!(matrix, values)
        end
    end
    return matrix
end


function save_temp(states, domain, iter, numeric_z_range, times)
    dir = "$domain" * "_temps"
    if !isdir(dir)
        mkdir(dir)
    end
    data = reduce(vcat, states)'
    data = reshape(data, length(numeric_z_range), length(times))
    file_path = joinpath(dir, "iter_$iter.csv")
    df = DataFrame(data, Symbol.(times))
    df = insertcols!(df, 1, :RowHeader => numeric_z_range)
    CSV.write(file_path, df)
end


function coupled_heat_equations()
    parameters = (
        h_atm=Float64(200),   # depth [m]
        h_oce=Float64(50),    # depth [m]
        n_atm=200,
        n_oce=50,
        k_atm=Float64(0.02364),
        k_oce=Float64(0.01),
        c_atm=Float64(1000),  # specific heat [J / kg / K]
        c_oce=Float64(4180),   # specific heat [J / kg / K]
        ρ_atm=Float64(1),     # density [kg / m3]
        ρ_oce=Float64(1000),  # density [kg / m3]
        u_atm=Float64(30),  # [m/s]
        u_oce=Float64(5),   #[m/s]
        C_AO=Float64(1e-5),
        C_AI=Float64(1e-5),
        C_OI=Float64(5e-5),
        T_atm_ini=Float64(280),   # initial temperature [K]
        T_oce_ini=Float64(288),   # initial temperature [K]
        T_ice_ini=Float64(272),       # temperature [K]
        a_i=Float64(0),           # ice area fraction [0-1]
        Δt_min=Float64(1.0),
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
        Δt_min=Float64(10.0),
        timerange=(Float64(0.0), Float64(1000)),
        Δt_coupler=Float64(1000.0),
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

    mesh_time_oce = CC.Meshes.IntervalMesh(time_points_oce_domain, nelems=Int(stepping.timerange[2] / stepping.Δt_min + 1))
    mesh_time_atm = CC.Meshes.IntervalMesh(time_points_atm_domain, nelems=Int(stepping.timerange[2] / stepping.Δt_min + 1))
    space_time_oce = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_time_oce)
    space_time_atm = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_time_atm)

    atmos_cache = (; parameters..., T_sfc=parameters.T_oce_ini .* CC.Fields.ones(space_time_oce), T_ice=T_ice_0)
    atmos_sim = atmos_init(stepping, T_atm_0, center_space_atm, atmos_cache)
    ocean_cache = (; parameters..., T_air=parameters.T_atm_ini .* CC.Fields.ones(space_time_atm), T_ice=T_ice_0)
    ocean_sim = ocean_init(stepping, T_oce_0, center_space_oce, ocean_cache)
    ice_cache = (; parameters...)
    ice_sim = ice_init(stepping, T_ice_0, point_space_ice, ice_cache)
    # println(parent(atmos_sim.integrator.p.T_sfc)[1])
    # println(typeof(atmos_sim.integrator.p))
    # println(fieldnames(typeof(atmos_sim.integrator.p)))

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

    solve_coupler!(cs, 1e-10, false, true)

end;

