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

function rename_file(cs::Interfacer.CoupledSimulation, iter, time, reverse=false)
    """When a file is saved, its always called the same thing.
    Had to rename it for each iteration to not overwrite"""

    for sim in cs.model_sims
        if !(Interfacer.name(sim) == "ConstantIce")
            original_file = joinpath(cs.dirs.artifacts, "checkpoint", "checkpoint_" * Interfacer.name(sim) * "_$time.hdf5")
            new_file = joinpath(cs.dirs.artifacts, "checkpoint", "checkpoint_" * Interfacer.name(sim) * "_$iter" * "_$time.hdf5")
            if !reverse
                mv(original_file, new_file, force=true)
            else
                mv(new_file, original_file, force=true)
            end
        end
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

function solve_coupler!(cs::Interfacer.CoupledSimulation, print_conv, plot_conv, return_conv)
    (; Δt_cpl, tspan) = cs

    cs.dates.date[1] = TimeManager.current_date(cs, tspan[begin])

    @info("Starting coupling loop")

    for t in ((tspan[begin]+Δt_cpl):Δt_cpl:tspan[end])

        time = Int(t - Δt_cpl)
        @info(cs.dates.date[1])
        iter = 1

        Checkpointer.checkpoint_sims(cs) # I had to remove nothing here
        rename_file(cs, 0, time)

        times = []
        numeric_z_range_ocean = []
        numeric_z_range_atmos = []

        # Save the z- and t-ranges
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
        bound_atmos_vals = Nothing
        bound_ocean_vals = Nothing

        # Should handle when it doesnt converge here (for example when delta z too small.)
        while true
            @info("Current iter: $(iter)")
            if iter > 10
                return NaN, NaN
            end
            # Update models
            Interfacer.step!(cs.model_sims.ice_sim, t)
            Interfacer.step!(cs.model_sims.ocean_sim, t)

            # Update atmosphere simulation
            ice_T = get_field(cs.model_sims.ice_sim, Val(:T_ice))
            ocean_states = copy(cs.model_sims.ocean_sim.integrator.sol.u)
            ocean_T = mean([ocean_state[end] for ocean_state in ocean_states])
            update_field!(cs.model_sims.atmos_sim, Val(:T_oce_sfc), ocean_T, Val(:T_ice), ice_T)

            # Step with atmosphere simulation and save atmosphere states
            Interfacer.step!(cs.model_sims.atmos_sim, t)
            atmos_states = copy(cs.model_sims.atmos_sim.integrator.sol.u)
            atmos_T = mean([atmos_state[1] for atmos_state in atmos_states])

            # Temperature values for this iteration.
            pre_bound_atmos_vals = bound_atmos_vals
            pre_bound_ocean_vals = bound_ocean_vals
            atmos_vals = extract_matrix(atmos_states, "atm")
            ocean_vals = extract_matrix(ocean_states, "oce")
            bound_atmos_vals = [atmos_val[1] for atmos_val in atmos_vals]
            bound_ocean_vals = [ocean_val[end] for ocean_val in ocean_vals]
            if iter > 1
                bound_errors_atm_iter = abs.(bound_atmos_vals .- pre_bound_atmos_vals)
                bound_errors_oce_iter = abs.(bound_ocean_vals .- pre_bound_ocean_vals)
                tols_atm = 100 * eps.(max.(abs.(bound_atmos_vals), abs.(pre_bound_atmos_vals)))
                tols_oce = 100 * eps.(max.(abs.(bound_ocean_vals), abs.(pre_bound_ocean_vals)))
                if all(bound_errors_atm_iter .< tols_atm)
                    println("Stopped at iter $iter for the atmosphere")
                    break
                elseif all(bound_errors_oce_iter .< tols_oce)
                    println("Stopped at iter $iter for the ocean")
                    break
                end
            end
            push!(atmos_vals_list, bound_atmos_vals)
            push!(ocean_vals_list, bound_ocean_vals)

            Checkpointer.checkpoint_sims(cs) # I had to remove nothing here

            rename_file(cs, iter, time)

            iter += 1
            restart_sims!(cs)
            reset_time!(cs, t - Δt_cpl)
            update_field!(cs.model_sims.ocean_sim, Val(:T_atm_sfc), atmos_T, Val(:T_ice), ice_T)
        end
        cs.dates.date[1] = TimeManager.current_date(cs, t)

        if print_conv || plot_conv || return_conv
            conv_fac_atm = []
            conv_fac_oce = []
            bound_error_atm = 0
            bound_error_oce = 0
            for i = 1:iter
                pre_bound_error_atm = bound_error_atm
                pre_bound_error_oce = bound_error_oce
                bound_error_atm = abs.(atmos_vals_list[i] .- atmos_vals_list[end])
                bound_error_oce = abs.(ocean_vals_list[i] .- ocean_vals_list[end])
                tols_atm = 100 * eps.(max.(abs.(atmos_vals_list[i]), abs.(atmos_vals_list[end])))
                tols_oce = 100 * eps.(max.(abs.(ocean_vals_list[i]), abs.(ocean_vals_list[end])))
                if all(isnan.(bound_error_atm) .|| (bound_error_atm .< tols_atm)) || all(isnan.(bound_error_oce) .|| (bound_error_oce .< tols_oce))
                    println("Too small values, cannot compute convergence factors.")
                    break
                end
                if i > 1
                    indices_atm = findall((pre_bound_error_atm .>= tols_atm) .& (bound_error_atm .>= tols_atm))
                    indices_oce = findall((pre_bound_error_oce .>= tols_oce) .& (pre_bound_error_oce .>= tols_oce))# This is maybe not good. Dont want the nominator to be less than the precision either, right. Should change and redo. Obs that the velocities are changed to -5:5!
                    conv_fac_atm_value = maximum(bound_error_atm[indices_atm] ./ pre_bound_error_atm[indices_atm])
                    conv_fac_oce_value = maximum(bound_error_oce[indices_oce] ./ pre_bound_error_oce[indices_oce])
                    push!(conv_fac_atm, conv_fac_atm_value)
                    push!(conv_fac_oce, conv_fac_oce_value)
                end
            end
            if print_conv
                println("Convergence factor atmosphere: $conv_fac_atm")
                println("Convergence factor ocean: $conv_fac_oce")
            end
            if plot_conv
                k_atm = 2:length(conv_fac_atm)+1
                k_oce = 2:length(conv_fac_oce)+1
                gr()
                scatter(k_atm, conv_fac_atm, label="atm", legend=:topright, color=:blue, markersize=5, xlabel="k", ylabel="ρₖ", ylim=(0, 0.0007))
                scatter!(k_oce, conv_fac_oce, label="oce", color=:green, markersize=5)
                display(current())

            end # Allow for computation, plot and print of convergence factor in this script.
            if return_conv
                return conv_fac_atm, conv_fac_oce
            end
        end
    end
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

function psi(xi, atm, stable, heat)
    if atm && stable && heat
        return -2 / 3 * (xi - (5 / 0.35)) * exp(-0.35 * xi) - (1 + (2 * xi / 3))^1.5 - (10 / 1.05) + 1
    elseif atm && stable && !heat
        return -2 / 3 * (xi - (5 / 0.35)) * exp(-0.35 * xi) - xi - (10 / 1.05)
    elseif atm && !stable && heat
        return 2 * log((1 + (1 - 16 * xi)^(1 / 2)) / 2)
    elseif atm && !stable && !heat
        x = (1 - 16 * xi)^(1 / 4)
        return pi / 2 - 2 * atan(x) + log((1 + x)^2 * (1 + x^2) / 8)
    elseif !atm && stable
        return 1 + 5 * xi
    elseif !atm && !stable
        if heat
            x = (1 - 25 * xi)^(1 / 3)
        else
            x = (1 - 14 * xi)^(1 / 3)
        end
        return sqrt(3) * (atan(sqrt(3)) - atan(1 / sqrt(3)) * (2 * x + 1)) + (3 / 2) * log((x^2 + x + 1) / 3)
    end
end

function coupled_heat_equations()
    # Parameters
    ρ_atm = Float64(1)
    ρ_oce = Float64(1000)
    c_atm = Float64(1000)
    c_oce = Float64(4180)
    lambda_u = sqrt(ρ_atm / ρ_oce)
    lambda_T = sqrt(ρ_atm / ρ_oce) * c_atm / c_oce
    nu_O = Float64(1e-6)
    nu_A = Float64(1.5 * 1e-5)
    mu = nu_O / nu_A
    kappa = Float64(0.4)

    # roughness lengths and C_OI
    z_0numA = Float64(10) # numerical zero height atmosphere [m]
    z_0numO = Float64(1) # numerical zero height atmosphere [m]
    z_ruAO = Float64(2 * 10^-4)
    z_rTAO = Float64(2 * 10^-4)
    z_rTAI = Float64(10^-3)
    C_OI = 5 * 1e-3


    # Variables
    a_i = 0.5
    # a_is = 0:0.01:1
    # a_is = 0.5
    # a_is = [0.1, 0.3, 0.5]

    t_max = 1000
    # t_maxs = [10, 100, 1000, 10000, 100000, 1000000]

    H_A = 210
    # H_As = [15, 50, 100, 200, 300, 500, 700, 1500] # Have to change n_atm here also
    H_O = 51
    # H_Os = [5, 10, 50, 100, 125, 150, 175, 200] # Have to change n_oce here also
    # n_atm = 200
    # n_atms = [20, 200, 2000, 20000]
    n_oce = 50
    # n_oces = [5, 27, 50, 275, 500, 2750, 5000, 27500, 50000]
    delta_z = 0.1:0.1:10
    n_atms = Int.(round.((H_A - z_0numA) ./ reverse(delta_z)))
    n_atms_unique = unique(n_atms)
    # delta_z = 0.001:0.001:0.1
    # n_oces = Int.(round.((H_O - z_0numO) ./ reverse(delta_z)))
    # n_oces_unique = unique(n_oces)
    # n_oces = Int.(round.(exp10.(range(log10(5000), log10(50000), length=10))))
    # n_atms = Int.(round.(exp10.(range(log10(40), log10(1000000), length=100))))

    L_AO = 50
    L_AI = 50
    L_OA = 100
    # L_AIs = [10, 50, 100, 125, 150, 175, 200]
    # C_AOs = [0.0007, 0.0008, 0.0009, 0.001, 0.0011, 0.0012, 0.0013]

    u_A = 5
    u_O = 1
    # u_As = [0, 0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
    # u_Os = [0, 0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]

    T_atm_ini = 267
    T_oce_ini = 276
    T_ice = 270
    # T_atm_inis = [260, 262, 264, 266, 268, 270, 273]
    # T_oce_inis = [270, 271, 272, 273, 274, 275, 276]
    # T_ices = [260, 262, 264, 266, 268, 270, 273]

    Δt_min = 100
    # n_t_oce = 1
    # n_t_atm = 500
    # Δt_min = 100
    n_t_oce = 10
    # n_t_atms = 1
    # n_t_atms = [1, 10, 100, 1000]
    # n_t_oces = [1, 10, 100, 1000]
    delta_t = 1:1:100
    n_t_atms = Int.(round.(Δt_min ./ reverse(delta_t)))
    n_t_atms_unique = unique(n_t_atms)
    # delta_t = 1:1:100
    # n_t_oces = Int.(round.(Δt_min ./ reverse(delta_t)))
    # n_t_oces_unique = unique(n_t_oces)
    # n_t_oces = Int.(round.(exp10.(range(log10(1), log10(10), length=10))))

    # conv_facs_atm = []
    # conv_facs_oce = []

    # conv_facs_atm = zeros(length(a_is), length(t_maxs))
    # conv_facs_oce = zeros(length(a_is), length(t_maxs))

    # conv_facs_atm = zeros(length(a_is), length(H_As))
    # conv_facs_oce = zeros(length(a_is), length(H_As))

    # conv_facs_atm = zeros(length(a_is), length(H_Os))
    # conv_facs_oce = zeros(length(a_is), length(H_Os))

    # conv_facs_atm = zeros(length(a_is), length(n_atms))
    # conv_facs_oce = zeros(length(a_is), length(n_atms))

    # conv_facs_atm = zeros(length(a_is), length(n_oces))
    # conv_facs_oce = zeros(length(a_is), length(n_oces))

    # conv_facs_atm = zeros(length(a_is), length(L_AIs))
    # conv_facs_oce = zeros(length(a_is), length(L_AIs))
    # C_AIs = zeros(length(a_is), length(L_AIs))

    # conv_facs_atm = zeros(length(a_is), length(C_AOs))
    # conv_facs_oce = zeros(length(a_is), length(C_AOs))

    # conv_facs_atm = zeros(length(a_is), length(u_As))
    # conv_facs_oce = zeros(length(a_is), length(u_As))

    # conv_facs_atm = zeros(length(a_is), length(u_Os))
    # conv_facs_oce = zeros(length(a_is), length(u_Os))

    # conv_facs_atm = zeros(length(a_is), length(T_atm_inis))
    # conv_facs_oce = zeros(length(a_is), length(T_atm_inis))

    # conv_facs_atm = zeros(length(a_is), length(T_oce_inis))
    # conv_facs_oce = zeros(length(a_is), length(T_oce_inis))

    # conv_facs_atm = zeros(length(a_is), length(T_ices))
    # conv_facs_oce = zeros(length(a_is), length(T_ices))

    # conv_facs_atm = zeros(length(a_is), length(n_t_atms))
    # conv_facs_oce = zeros(length(a_is), length(n_t_atms))

    # conv_facs_atm = zeros(length(a_is), length(n_t_oces))
    # conv_facs_oce = zeros(length(a_is), length(n_t_oces))

    conv_facs_atm = zeros(length(n_t_atms_unique), length(n_atms_unique))
    conv_facs_oce = zeros(length(n_t_atms_unique), length(n_atms_unique))

    # conv_facs_atm = zeros(length(n_t_oces_unique), length(n_oces_unique))
    # conv_facs_oce = zeros(length(n_t_oces_unique), length(n_oces_unique))
    # for (k, a_i) in enumerate(a_is)
    for (k, n_t_atm) in enumerate(n_t_atms_unique)
        for (j, n_atm) in enumerate(n_atms_unique) # remove if you only want rho as function of a_i
            println("n_t_oce=$n_t_oce")
            # n_atm = H_A - 10 #If varing H_A to not vary deltaz
            # n_oce = H_O - 1 #If varing H_A to not vary deltaz
            # Variables
            a_i = Float64(a_i)
            z_ruAI = Float64(maximum([10^-3, 0.93 * 10^-3 * (1 - a_i) + 6.05 * 10^-3 * exp(-17 * (a_i - 0.5)^2)]))
            L_AO = Float64(L_AO)
            L_AI = Float64(L_AI)
            L_OA = Float64(L_OA)
            if L_AO > 0
                stable_atm = true
            else
                stable_atm = false
            end
            if L_OA > 0
                stable_oce = true
            else
                stable_oce = false
            end
            C_AO = kappa^2 / ((log(z_0numA / z_ruAO) - psi(z_0numA / L_AO, true, stable_atm, false) + lambda_u * (log(lambda_u * z_0numO / (z_ruAO * mu)) - psi(z_0numO / L_OA, false, stable_oce, false))) * (log(z_0numA / z_rTAO) - psi(z_0numA / L_AO, true, stable_atm, true) + lambda_T * (log(lambda_T * z_0numO / (z_rTAO * mu)) - psi(z_0numO / L_OA, false, stable_oce, true))))
            C_AI = kappa^2 / (log(z_0numA / z_ruAI) - psi(z_0numA / L_AI, true, stable_atm, false)) * (log(z_0numA / z_rTAI) - psi(z_0numA / L_AI, true, stable_atm, true))
            # C_AIs[k, j] = C_AI
            parameters = (
                h_atm=Float64(H_A),   # depth [m]
                h_oce=Float64(H_O),    # depth [m]
                n_atm=n_atm,
                n_oce=n_oce,
                k_atm=Float64(0.02364),
                k_oce=Float64(0.01),
                c_atm=c_atm,  # specific heat [J / kg / K]
                c_oce=c_oce,   # specific heat [J / kg / K]
                ρ_atm=ρ_atm,     # density [kg / m3]
                ρ_oce=ρ_oce,  # density [kg / m3]
                u_atm=Float64(u_A),  # [m/s]
                u_oce=Float64(u_O),   #[m/s]
                C_AO=Float64(C_AO),
                C_AI=Float64(C_AI),
                C_OI=Float64(C_OI),
                T_atm_ini=Float64(T_atm_ini),   # initial temperature [K]
                T_oce_ini=Float64(T_oce_ini),   # initial temperature [K]
                T_ice_ini=Float64(T_ice),       # temperature [K]
                a_i=a_i,           # ice area fraction [0-1]
            )

            context = CC.ClimaComms.context()
            device = CC.ClimaComms.device(context)

            # initialize models
            domain_atm = CC.Domains.IntervalDomain(
                CC.Geometry.ZPoint{Float64}(z_0numA),
                CC.Geometry.ZPoint{Float64}(parameters.h_atm);
                boundary_names=(:bottom, :top),
            )
            mesh_atm = CC.Meshes.IntervalMesh(domain_atm, nelems=parameters.n_atm)
            center_space_atm = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_atm)

            domain_oce = CC.Domains.IntervalDomain(
                CC.Geometry.ZPoint{Float64}(-parameters.h_oce),
                CC.Geometry.ZPoint{Float64}(-z_0numO);
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
                Δt_min=Float64(Δt_min),
                timerange=(Float64(0.0), Float64(t_max)),
                Δt_coupler=Float64(t_max),
                odesolver=CTS.ExplicitAlgorithm(CTS.RK4()),
                nsteps_atm=n_t_atm,
                nsteps_oce=n_t_oce,
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

            conv_fac_atm, conv_fac_oce = solve_coupler!(cs, true, false, true)
            # solve_coupler!(cs, false, true, false)
            if !(conv_fac_atm isa AbstractArray) && isnan(conv_fac_atm)
                conv_facs_atm[k, j] = Inf
            elseif !isempty(conv_fac_atm)
                # push!(conv_facs_atm, conv_fac_atm[1])
                conv_facs_atm[k, j] = conv_fac_atm[1]
            else
                # push!(conv_facs_atm, NaN)
                conv_facs_atm[k, j] = NaN
            end
            if !(conv_fac_oce isa AbstractArray) && isnan(conv_fac_oce)
                conv_facs_oce[k, j] = Inf
            elseif !isempty(conv_fac_oce)
                # push!(conv_facs_oce, conv_fac_oce[1])
                conv_facs_oce[k, j] = conv_fac_oce[1]
            else
                # push!(conv_facs_oce, NaN)
                conv_facs_oce[k, j] = NaN
            end
        end
    end

    mapping_vals_atm = Dict((x, y) => conv_facs_atm[findfirst(==(x), n_t_atms_unique), findfirst(==(y), n_atms_unique)]
                            for x in n_t_atms, y in n_atms)

    mapping_vals_oce = Dict((x, y) => conv_facs_oce[findfirst(==(x), n_t_atms_unique), findfirst(==(y), n_atms_unique)]
                            for x in n_t_atms, y in n_atms)

    println(mapping_vals_oce)

    conv_facs_atm = reshape([mapping_vals_atm[(x, y)] for x in n_t_atms, y in n_atms], (length(n_t_atms), length(n_atms)))
    conv_facs_oce = reshape([mapping_vals_oce[(x, y)] for x in n_t_atms, y in n_atms], (length(n_t_atms), length(n_atms)))

    # Print conv factors
    println(size(conv_facs_atm))
    println(conv_facs_oce)

    # Set appropriate backend for plotting
    # plotly()
    gr()

    # Plot conv factor as a function of a_i
    # plot(a_is[1:length(conv_facs_atm)], log.(conv_facs_atm), label="atm", color=:blue, markersize=5, xlabel="aᴵ", ylabel="log(ρ)", legend=:topright, ylim=(-20, -7))
    # scatter!(a_is[1:length(conv_facs_atm)], log.(conv_facs_atm), label="", color=:blue, markersize=5, xlabel="aᴵ", ylabel="log(ρ)", legend=:topright)
    # plot!(a_is[1:length(conv_facs_oce)], log.(conv_facs_oce), label="oce", color=:green, markersize=5, xlabel="aᴵ", ylabel="log(ρ)", legend=:topright)
    # scatter!(a_is[1:length(conv_facs_oce)], log.(conv_facs_oce), label="", color=:green, markersize=5, xlabel="aᴵ", ylabel="log(ρ)", legend=:topright)


    # plot conv factor for different a_i and wrt a parameter. Some special treatment
    # for divergence has been added if the variable is delta t or delta z
    # colors = [:blue, :red, :green]
    # for (i, ai) in enumerate(a_is)
    #     # x = reverse((H_O - z_0numO) ./ n_oces) # For delta_z
    #     # x = reverse(Δt_min ./ n_t_oces) # For delta_t
    #     y = reverse(conv_facs_atm[i, :])
    #     # y = reverse(conv_facs_oce[i, :])
    #     x_numeric = [x[j] for j in 1:length(y) if !isinf(y[j])]
    #     y_numeric = [y[j] for j in 1:length(y) if !isinf(y[j])]
    #     x_text = [x[j] for j in 1:length(y) if isinf(y[j])]
    #     y_text = [5e-5 + (i - 2) * 6e-6 for _ in x_text] # For delta_z_A
    #     # y_text = [5e-6 + (i - 2) * 6e-6 for _ in x_text] # For delta_z_O
    #     # x = reverse((H_O - z_0numO) ./ n_oces)
    #     x_axis_label = "Δzᴼ"
    #     if i == 1
    #         # scatter(x, conv_facs_atm[i, :], xticks=x, xscale=:log10, yformatter=:scientific, label="aᴵ=$ai", markersize=5, xlabel=x_axis_label, ylabel="ρᴬ", legend=:topright, color=colors[i])
    #         # plot!(x, conv_facs_atm[i, :], label="", linewidth=2, color=colors[i])
    #         # scatter(x, conv_facs_oce[i, :], xticks=x, xscale=:log10, yformatter=:scientific, label="aᴵ=$ai", markersize=5, xlabel=x_axis_label, ylabel="ρᴼ", legend=:topright, color=colors[i])
    #         # plot!(x, conv_facs_oce[i, :], label="", linewidth=2, color=colors[i])

    #         scatter(x_numeric, y_numeric, xticks=x, xscale=:log10, yformatter=:scientific, label="aᴵ=$ai", markersize=5, xlabel=x_axis_label, ylabel="ρᴬ", legend=:topright, color=colors[i])
    #         plot!(x_numeric, y_numeric, label="", linewidth=2, color=colors[i])
    #         for (k, txt) in enumerate(y_text)
    #             annotate!(x_text[k], txt, text("Divergence", 12, colors[i], :center, rotation=90))
    #         end
    #     else
    #         # scatter!(x, conv_facs_atm[i, :], label="aᴵ=$ai", markersize=5, color=colors[i])
    #         # plot!(x, conv_facs_atm[i, :], label="", linewidth=2, color=colors[i])
    #         # scatter!(x, conv_facs_oce[i, :], label="aᴵ=$ai", markersize=5, color=colors[i])
    #         # plot!(x, conv_facs_oce[i, :], label="", linewidth=2, color=colors[i])
    #         scatter!(x_numeric, y_numeric, xticks=x, xscale=:log10, yformatter=:scientific, label="aᴵ=$ai", markersize=5, xlabel=x_axis_label, ylabel="ρᴬ", legend=:topright, color=colors[i])
    #         plot!(x_numeric, y_numeric, label="", linewidth=2, color=colors[i])
    #         for (k, txt) in enumerate(y_text)
    #             println(x_text[k])
    #             println(txt)
    #             annotate!(x_text[k], txt, text("Divergence", 12, colors[i], :center, rotation=90))
    #         end
    #     end
    # end


    # # Plot convergence factor as a function of deltaz delta t quotient
    # Δzᴼs = reverse((H_O - z_0numO) ./ n_oces) # For delta_z
    # quotient = Δzᴼs ./ (Δt_min / n_t_oce)
    # conv_facs_atm[.!isinf.(conv_facs_atm)] .= 0
    # conv_facs_atm[isinf.(conv_facs_atm)] .= 1
    # heatmap(conv_facs_atm', xticks=(1:20:length(a_is), round.(a_is[1:20:length(a_is)], digits=1)), yticks=(1:20:length(quotient), round.(quotient[1:20:length(quotient)], digits=4)), color=:viridis,
    #     xlabel="aᴵ", ylabel="Δzᴼ/Δtᴼ",
    #     clims=(-1, 1))  # Set color limits    println(conv_facs_atm)
    # # plot3d(a_is, quotient, conv_facs_atm', label="atm", color=:blue, markersize=5, xlabel="aᵢ", ylabel="", zlabel="ρ", legend=:topright)
    # # # surface!(a_is[1:length(conv_facs_oce[:, 1])], t_maxs[1:length(conv_facs_oce[1, :])], conv_facs_oce', label="oce", color=:green, markersize=5)
    # display(current())

    # Plot convergence factor as a function of deltaz and delta t
    # Δzᴼs = reverse((H_O - z_0numO) ./ n_oces) # For delta_z
    # Δtᴼs = reverse(Δt_min ./ n_t_oces) # For delta_t

    conv_facs_atm[.!isinf.(conv_facs_atm)] .= 0
    conv_facs_atm[isinf.(conv_facs_atm)] .= 1
    conv_facs_atm = reverse(conv_facs_atm, dims=(1, 2))
    heatmap(conv_facs_atm', xticks=(1:20:length(delta_t), delta_t[1:20:length(delta_t)]), yticks=(1:20:length(delta_z), round.(delta_z[1:20:length(delta_z)], digits=3)), color=:viridis,
        xlabel="Δtᴬ", ylabel="Δzᴬ",
        clims=(0, 1))  # Set color limits    println(conv_facs_atm)
    # plot3d(a_is, quotient, conv_facs_atm', label="atm", color=:blue, markersize=5, xlabel="aᵢ", ylabel="", zlabel="ρ", legend=:topright)
    # # surface!(a_is[1:length(conv_facs_oce[:, 1])], t_maxs[1:length(conv_facs_oce[1, :])], conv_facs_oce', label="oce", color=:green, markersize=5)
    println(n_t_atms)
    println(n_atms)
    # println(Δtᴼs)
    # println(Δzᴼs)
    display(current())


end;

