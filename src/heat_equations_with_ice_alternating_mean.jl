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
        stopped_at_nan = false
        starting_temp_oce = cs.model_sims.atmos_sim.params.T_sfc[]
        starting_temp_atm = cs.model_sims.ocean_sim.params.T_air[]
        starting_temp_ice = cs.model_sims.atmos_sim.params.T_ice[]
        upper_limit_temp = maximum([starting_temp_oce, starting_temp_atm, starting_temp_ice])
        lower_limit_temp = minimum([starting_temp_oce, starting_temp_atm, starting_temp_ice])

        # Should handle when it doesnt converge here (for example when delta z too small.)
        while true
            @info("Current iter: $(iter)")
            if iter == 10
                println("Stopped at iter 10 due too slow convergence or coupling divergence")
                break
            end
            # Update models
            Interfacer.step!(cs.model_sims.ice_sim, t)
            Interfacer.step!(cs.model_sims.ocean_sim, t)

            # Update atmosphere simulation
            ice_T = get_field(cs.model_sims.ice_sim, Val(:T_ice))
            ocean_states = copy(cs.model_sims.ocean_sim.integrator.sol.u)
            ocean_T = mean([ocean_state[end] for ocean_state in ocean_states])
            println(ocean_T)
            println(ice_T)
            update_field!(cs.model_sims.atmos_sim, Val(:T_oce_sfc), ocean_T, Val(:T_ice), ice_T)

            # Step with atmosphere simulation and save atmosphere states
            Interfacer.step!(cs.model_sims.atmos_sim, t)
            atmos_states = copy(cs.model_sims.atmos_sim.integrator.sol.u)
            atmos_T = mean([atmos_state[1] for atmos_state in atmos_states])
            println(atmos_T)

            # Temperature values for this iteration.
            pre_bound_atmos_vals = bound_atmos_vals
            pre_bound_ocean_vals = bound_ocean_vals
            atmos_vals = extract_matrix(atmos_states, "atm")
            ocean_vals = extract_matrix(ocean_states, "oce")
            # bound_atmos_vals = [atmos_val[1] for atmos_val in atmos_vals]
            # bound_ocean_vals = [ocean_val[end] for ocean_val in ocean_vals]
            bound_atmos_vals = atmos_vals[1, :]
            bound_ocean_vals = ocean_vals[end, :]
            if (any(isnan, atmos_vals) || any(isnan, ocean_vals) || maximum(ocean_vals) > upper_limit_temp || minimum(ocean_vals) < lower_limit_temp || maximum(atmos_vals) > upper_limit_temp || minimum(atmos_vals) < lower_limit_temp) && iter == 1
                println("stopped due to unstable model")
                stopped_at_nan = true
                break
            elseif (any(isnan, bound_atmos_vals) || any(isnan, bound_ocean_vals)) && iter > 1
                println("stopped due to coupling divergence")
                break
            end
            if iter > 1
                bound_errors_atm_iter = abs.(bound_atmos_vals .- pre_bound_atmos_vals)
                bound_errors_oce_iter = abs.(bound_ocean_vals .- pre_bound_ocean_vals)
                tols_atm = 100 * eps.(max.(abs.(bound_atmos_vals), abs.(pre_bound_atmos_vals)))
                tols_oce = 100 * eps.(max.(abs.(bound_ocean_vals), abs.(pre_bound_ocean_vals)))
                if all(bound_errors_atm_iter .< tols_atm) || all(bound_errors_oce_iter .< tols_oce)
                    println("stopped at iter $iter")
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

        if stopped_at_nan
            return NaN, NaN
        else
            end_of_loop = length(atmos_vals_list) - 1
        end
        if print_conv || plot_conv || return_conv
            conv_fac_atm = []
            conv_fac_oce = []
            bound_error_atm = 0
            bound_error_oce = 0
            for i = 1:end_of_loop
                pre_bound_error_atm = bound_error_atm
                pre_bound_error_oce = bound_error_oce
                bound_error_atm = abs.(atmos_vals_list[i] .- atmos_vals_list[end])
                bound_error_oce = abs.(ocean_vals_list[i] .- ocean_vals_list[end])
                tols_atm = 100 * eps.(max.(abs.(atmos_vals_list[i]), abs.(atmos_vals_list[end])))
                tols_oce = 100 * eps.(max.(abs.(ocean_vals_list[i]), abs.(ocean_vals_list[end])))
                if i > 1
                    indices_atm = findall((pre_bound_error_atm .>= tols_atm) .& (bound_error_atm .>= tols_atm))
                    indices_oce = findall((pre_bound_error_oce .>= tols_oce) .& (pre_bound_error_oce .>= tols_oce))# This is maybe not good. Dont want the nominator to be less than the precision either, right. Should change and redo. Obs that the velocities are changed to -5:5!
                    conv_fac_atm_value = maximum(bound_error_atm[indices_atm] ./ pre_bound_error_atm[indices_atm]) # Change to sum(.)/sum(.)
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
    return hcat(matrix...)
end

function psi(xi, atm, stable, heat)
    if stable
        if atm && heat
            return -2 / 3 * (xi - (5 / 0.35)) * exp(-0.35 * xi) - (1 + (2 * xi / 3))^1.5 - (10 / 1.05) + 1
        elseif atm && !heat
            return -2 / 3 * (xi - (5 / 0.35)) * exp(-0.35 * xi) - xi - (10 / 1.05)
        elseif !atm
            return -5 * xi
        end
    else
        if !atm && heat
            x = (1 - 25 * xi)^(1 / 3)
            return sqrt(3) * (atan(np.sqrt(3)) - atan(1 / np.sqrt(3) * (2 * x + 1))) + (3 / 2) * log((x^2 + x + 1) / 3)
        elseif !atm && !heat
            x = (1 - 14 * xi)^(1 / 3)
            return sqrt(3) * (atan(np.sqrt(3)) - atan(1 / np.sqrt(3) * (2 * x + 1))) + (3 / 2) * log((x^2 + x + 1) / 3)
        elseif atm && heat
            return 2 * np.log((1 + (1 - 16 * xi)^(1 / 2)) / 2)
        else
            x = (1 - 16 * xi)^(1 / 4)
            return np.pi / 2 - 2 * np.atan(x) + np.log((1 + x)^2 * (1 + x^2) / 8)
        end
    end
end

function define_realistic_vals()

    physical_values = Dict(
        :a_i => Float64(0),
        :ρ_atm => Float64(1.225),
        :ρ_oce => Float64(1000.0),
        :c_atm => Float64(1005.0),
        :c_oce => Float64(4182.0),
        :lambda_u => Float64(sqrt(1.225 / 1000.0)),
        :lambda_T => Float64(sqrt(1.225 / 1000.0) * 1005.0 / 4182.0),
        :nu_O => Float64(1e-6),
        :nu_A => Float64(1.5e-5),
        :mu => Float64(1e-6 / 1.5e-5),
        :kappa => Float64(0.4),
        :k_atm => Float64(0.02364),
        :k_oce => Float64(0.58),
        :alpha_o => Float64(0.58 / (1000.0 * 4182.0)),
        :alpha_a => Float64(0.02364 / (1.225 * 1005.0)),
        :alpha_eos => Float64(1.8e-4),
        :z_0numA => Float64(10.0),
        :z_0numO => Float64(1.0),
        :z_ruAO => Float64(2e-4),
        :z_ruAI => nothing,
        :z_rTAO => Float64(2e-4),
        :z_rTAI => Float64(1e-3),
        :C_OI => Float64(5e-3),
        :h_oce => Float64(51.0),
        :h_atm => Float64(210.0),
        :u_atm => Float64(5.0),
        :u_oce => Float64(1.0),
        :L_AO => Float64(50.0),
        :L_AI => Float64(50.0),
        :L_OA => nothing,  # Placeholder since it depends on `T_atm_ini`
        :C_AO => nothing,  # Placeholder for computation
        :C_AI => nothing,  # Placeholder for computation
        :T_atm_ini => Float64(267.0),
        :T_oce_ini => Float64(276.0),
        :T_ice_ini => Float64(270.0),
        :t_max => Float64(1000.0),
        :Δt_min => Float64(10.0),
        :n_t_atm => 50,
        :n_t_oce => 1,
        :n_atm => 200,
        :n_oce => 50
    )

    physical_values[:L_OA] = physical_values[:lambda_u]^2 / (physical_values[:T_atm_ini] * physical_values[:alpha_eos] * physical_values[:lambda_T])
    physical_values[:C_AO] = physical_values[:kappa]^2 / ((log(physical_values[:z_0numA] / physical_values[:z_ruAO]) - psi(physical_values[:z_0numA] / physical_values[:L_AO], true, physical_values[:L_AO] > 0, false) + physical_values[:lambda_u] * (log(physical_values[:lambda_u] * physical_values[:z_0numO] / (physical_values[:z_ruAO] * physical_values[:mu])) - psi(physical_values[:z_0numO] / physical_values[:L_OA], false, physical_values[:L_OA] > 0, false))) * (log(physical_values[:z_0numA] / physical_values[:z_rTAO]) - psi(physical_values[:z_0numA] / physical_values[:L_AO], true, physical_values[:L_AO] > 0, true) + physical_values[:lambda_T] * (log(physical_values[:lambda_T] * physical_values[:z_0numO] / (physical_values[:z_rTAO] * physical_values[:mu])) - psi(physical_values[:z_0numO] / physical_values[:L_OA], false, physical_values[:L_OA] > 0, true))))
    physical_values[:w_min] = pi / physical_values[:t_max]
    z_ruAI = Float64(max(1e-3, 0.93e-3 * (1 - physical_values[:a_i]) + 6.05e-3 * exp(-17 * (physical_values[:a_i] - 0.5)^2)))  # Placeholder since it depends on `a_i`
    physical_values[:z_ruAI] = z_ruAI
    physical_values[:C_AI] = physical_values[:kappa]^2 / (log(physical_values[:z_0numA] / physical_values[:z_ruAI]) - psi(physical_values[:z_0numA] / physical_values[:L_AI], true, physical_values[:L_AI] > 0, false)) * (log(physical_values[:z_0numA] / physical_values[:z_rTAI]) - psi(physical_values[:z_0numA] / physical_values[:L_AI], true, physical_values[:L_AI] > 0, true))
    return physical_values
end

function update_physical_values(a_i, physical_values)
    physical_values[:a_i] = a_i
    z_ruAI = Float64(max(1e-3, 0.93e-3 * (1 - a_i) + 6.05e-3 * exp(-17 * (a_i - 0.5)^2)))  # Placeholder since it depends on `a_i`
    physical_values[:z_ruAI] = z_ruAI
    physical_values[:C_AI] = physical_values[:kappa]^2 / (log(physical_values[:z_0numA] / physical_values[:z_ruAI]) - psi(physical_values[:z_0numA] / physical_values[:L_AI], true, physical_values[:L_AI] > 0, false)) * (log(physical_values[:z_0numA] / physical_values[:z_rTAI]) - psi(physical_values[:z_0numA] / physical_values[:L_AI], true, physical_values[:L_AI] > 0, true))
    return physical_values
end

function compute_actual_rho(physical_values)
    sigma_o = sqrt(abs(physical_values[:w_min]) / (2 * physical_values[:alpha_o])) * (1 + im)
    sigma_a = sqrt(abs(physical_values[:w_min]) / (2 * physical_values[:alpha_a])) * (1 + im)
    eta_AO = physical_values[:C_AO] * abs(physical_values[:u_atm] - physical_values[:u_oce]) * physical_values[:ρ_atm] * physical_values[:c_atm]
    eta_OI = physical_values[:C_OI] * abs(physical_values[:u_oce]) * physical_values[:ρ_oce] * physical_values[:c_oce]
    eta_AI = physical_values[:C_AI] * abs(physical_values[:u_atm]) * physical_values[:ρ_atm] * physical_values[:c_atm]
    return abs((1 - physical_values[:a_i])^2 * eta_AO^2 / ((physical_values[:k_oce] * sigma_o * (1 / tanh(sigma_o * (physical_values[:h_oce] - physical_values[:z_0numO]))) + (1 - physical_values[:a_i]) * eta_AO + physical_values[:a_i] * eta_OI) * (physical_values[:k_atm] * sigma_a * (1 / tanh(sigma_a * (physical_values[:h_atm] - physical_values[:z_0numA]))) + (1 - physical_values[:a_i]) * eta_AO + physical_values[:a_i] * eta_AI)))
end

function get_coupled_sim(physical_values)
    physical_keys = [
        :h_atm, :h_oce, :n_atm, :n_oce, :k_atm, :k_oce, :c_atm, :c_oce,
        :ρ_atm, :ρ_oce, :u_atm, :u_oce, :C_AO, :C_AI, :C_OI, :T_atm_ini,
        :T_oce_ini, :T_ice_ini, :a_i
    ]
    parameters = NamedTuple(Dict(
        key => physical_values[key] for key in physical_keys
    ))

    context = CC.ClimaComms.context()
    device = CC.ClimaComms.device(context)

    # initialize models
    domain_atm = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(physical_values[:z_0numA]),
        CC.Geometry.ZPoint{Float64}(parameters.h_atm);
        boundary_names=(:bottom, :top),
    )
    mesh_atm = CC.Meshes.IntervalMesh(domain_atm, nelems=parameters.n_atm)
    center_space_atm = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_atm)

    domain_oce = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{Float64}(-parameters.h_oce),
        CC.Geometry.ZPoint{Float64}(-physical_values[:z_0numO]);
        boundary_names=(:bottom, :top),
    )
    mesh_oce = CC.Meshes.IntervalMesh(domain_oce, nelems=parameters.n_oce)
    center_space_oce = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_oce)

    # Adding height for the ice
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
        Δt_min=Float64(physical_values[:Δt_min]),
        timerange=(Float64(0.0), Float64(physical_values[:t_max])),
        Δt_coupler=Float64(physical_values[:t_max]),
        odesolver=CTS.ExplicitAlgorithm(CTS.RK4()),
        nsteps_atm=physical_values[:n_t_atm],
        nsteps_oce=physical_values[:n_t_oce],
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
    return cs
end

function compute_numerical_conv_fac(physical_values; return_conv_facs=true, plot_conv_facs=false, print_conv_facs=false)
    cs = get_coupled_sim(physical_values)
    if return_conv_facs
        conv_fac_atm, conv_fac_oce = solve_coupler!(cs, false, false, true)
        return conv_fac_atm, conv_fac_oce
    elseif plot_conv_facs
        solve_coupler!(cs, false, true, false)
    elseif print_conv_facs
        solve_coupler!(cs, true, false, false)
    end

end

function get_conv_facs_one_variable(physical_values, var1s, var1_name; a_i=0, analytic=false, log_scale=false)
    # Getting physical values
    if !(var1_name == "a_i")
        update_physical_values(a_i, physical_values)
    end

    # Introducing the convergence factors
    conv_facs_atm = zeros(length(var1s))
    conv_facs_oce = zeros(length(var1s))
    if analytic
        if log_scale
            variable1_range = exp10.(range(log10(var1s[1]), log10(var1s[end]), length=100))
        else
            variable1_range = range(var1s[1], var1s[end], length=100)
        end
        conv_facs_analytic = zeros(length(variable1_range))
    end

    # Looping over the variables
    for (k, var1) in enumerate(var1s)
        # Update physical values based on variable 1
        if var1_name == "a_i"
            update_physical_values(var1, physical_values)
        else
            physical_values[Symbol(var1_name)] = var1
        end
        # n_atm = H_A - 10 #If varing H_A to not vary deltaz
        # n_oce = H_O - 1 #If varing H_A to not vary deltaz
        # C_AIs[k, j] = C_AI
        conv_fac_atm, conv_fac_oce = compute_numerical_conv_fac(physical_values)
        if !(conv_fac_atm isa AbstractArray) && isnan(conv_fac_atm)
            conv_facs_atm[k] = Inf
        elseif !isempty(conv_fac_atm)
            conv_facs_atm[k] = conv_fac_atm[1]
        else
            conv_facs_atm[k] = NaN
        end
        if !(conv_fac_oce isa AbstractArray) && isnan(conv_fac_oce)
            conv_facs_oce[k] = Inf
        elseif !isempty(conv_fac_oce)
            conv_facs_oce[k] = conv_fac_oce[1]
        else
            conv_facs_oce[k] = NaN
        end
    end
    if analytic
        for (k, var1) in enumerate(variable1_range)
            if var1_name == "t_max"
                physical_values[:w_min] = pi / var1
            end
            physical_values[Symbol(var1_name)] = var1
            conv_fac_analytic = compute_actual_rho(physical_values)
            conv_facs_analytic[k] = conv_fac_analytic
        end
    end
    if analytic
        return conv_facs_atm, conv_facs_oce, variable1_range, conv_facs_analytic
    else
        return conv_facs_atm, conv_facs_oce
    end
end

function get_conv_facs_two_variables(physical_values, var1s, var2s, var1_name, var2_name; a_i=0, analytic=false, log_scale=false)
    """ 
    If a_i is sent in as a parameter, it should be the first one, 
    only allows for deltat and deltaz apart from having a_i and another variable.
    """

    # Getting physical values
    if !(var1_name == "a_i")
        # Updating physical values based on a_i if it is not a variable
        update_physical_values(a_i, physical_values)
    end

    # Introducing the convergence factors
    conv_facs_atm = zeros(length(var1s), length(var2s))
    conv_facs_oce = zeros(length(var1s), length(var2s))
    if analytic
        if log_scale
            variable2_range = exp10.(range(log10(var2s[1]), log10(var2s[end]), length=100))
        else
            variable2_range = range(var2s[1], var2s[end], length=100)
        end
        conv_facs_analytic = zeros(length(var1s), length(variable2_range))
    end

    if var2_name == "L_AI"
        C_AIs = zeros(length(var1s), length(var2s))
        C_AIs_analytic = zeros(length(var1s), length(variable2_range))
    end

    # Looping over the variables
    for (k, var1) in enumerate(var1s)
        # Update physical values based on variable 1
        if var1_name == "a_i"
            update_physical_values(var1, physical_values)
        else
            physical_values[Symbol(var1_name)] = var1
        end

        for (j, var2) in enumerate(var2s)
            physical_values[Symbol(var2_name)] = var2
            if var2_name == "L_AI"
                update_physical_values(var1, physical_values)
            end

            # Special treatments
            if var2_name == "h_atm" || var1_name == "h_atm"
                physical_values[:n_atm] = physical_values[:h_atm] - 10
            end
            if var2_name == "h_oce" || var1_name == "h_oce"
                physical_values[:n_oce] = physical_values[:h_oce] - 1
            end
            if var2_name == "L_AI"
                C_AIs[k, j] = physical_values[:C_AI]
            end

            # Compute convergence factor
            conv_fac_atm, conv_fac_oce = compute_numerical_conv_fac(physical_values)

            # Append convergence factor
            if !(conv_fac_atm isa AbstractArray) && isnan(conv_fac_atm)
                conv_facs_atm[k, j] = Inf
            elseif !isempty(conv_fac_atm)
                conv_facs_atm[k, j] = conv_fac_atm[1]
            else
                conv_facs_atm[k, j] = NaN
            end
            if !(conv_fac_oce isa AbstractArray) && isnan(conv_fac_oce)
                conv_facs_oce[k, j] = Inf
            elseif !isempty(conv_fac_oce)
                conv_facs_oce[k, j] = conv_fac_oce[1]
            else
                conv_facs_oce[k, j] = NaN
            end
        end

        # Compute analytical convergence factors for variable two. Only done in the case of a_i and another variable
        if analytic
            for (j, var2_analytic) in enumerate(variable2_range)
                if var2_name == "t_max"
                    physical_values[:w_min] = pi / var2_analytic
                end
                physical_values[Symbol(var2_name)] = var2_analytic
                if var2_name == "L_AI"
                    update_physical_values(var1, physical_values)
                    C_AIs_analytic[k, j] = physical_values[:C_AI]
                end
                conv_fac_analytic = compute_actual_rho(physical_values)
                conv_facs_analytic[k, j] = conv_fac_analytic
            end
        end
    end

    if !(var2_name == "L_AI")
        if analytic
            return conv_facs_atm, conv_facs_oce, variable2_range, conv_facs_analytic
        else
            return conv_facs_atm, conv_facs_oce
        end
    else
        if analytic
            return conv_facs_atm, conv_facs_oce, variable2_range, conv_facs_analytic, C_AIs, C_AIs_analytic
        else
            return conv_facs_atm, conv_facs_oce, C_AIs, C_AIs_analytic
        end
    end
end

function one_iter_update(physical_values, var1s, var2s, var1_name, var2_name)
    conv_facs_atm = zeros(length(var1s), length(var2s))
    conv_facs_oce = zeros(length(var1s), length(var2s))
    for (i, var1) in enumerate(var1s)
        for (j, var2) in enumerate(var2s)
            physical_values[Symbol(var1_name)] = var1
            physical_values[Symbol(var2_name)] = var2
            cs = get_coupled_sim(physical_values)
            starting_temp_oce = cs.model_sims.atmos_sim.params.T_sfc[]
            starting_temp_atm = cs.model_sims.ocean_sim.params.T_air[]
            starting_temp_ice = cs.model_sims.atmos_sim.params.T_ice[]
            upper_limit_temp = maximum([starting_temp_oce, starting_temp_atm, starting_temp_ice])
            lower_limit_temp = minimum([starting_temp_oce, starting_temp_atm, starting_temp_ice])
            Interfacer.step!(cs.model_sims.ocean_sim, physical_values[:t_max])
            Interfacer.step!(cs.model_sims.atmos_sim, physical_values[:t_max])
            ocean_states = copy(cs.model_sims.ocean_sim.integrator.sol.u)
            atmos_states = copy(cs.model_sims.atmos_sim.integrator.sol.u)

            ocean_vals = extract_matrix(ocean_states, "oce")
            atmos_vals = extract_matrix(atmos_states, "atm")
            if (any(isnan, atmos_vals) || maximum(atmos_vals) > upper_limit_temp || minimum(atmos_vals) < lower_limit_temp)
                conv_facs_atm[i, j] = Inf
            else
                conv_facs_atm[i, j] = NaN
            end
            if any(isnan, ocean_vals) || maximum(ocean_vals) > upper_limit_temp || minimum(ocean_vals) < lower_limit_temp
                conv_facs_oce[i, j] = Inf
            else
                conv_facs_oce[i, j] = NaN
            end

        end
    end
    return conv_facs_atm, conv_facs_oce
end
function coupled_heat_equations()
    # To run them separately: 
    # Create a method that creates the coupled simulation. Step only for one of the domains separately. 
    # See if the instability behaviour changes compared to when i had coupling. Do for different ice concentrations? Just do interfacer.step and check the same conditions as before.
    # Then i get NaNs if it is problematic. and i can check
    # Get physical values
    physical_values = define_realistic_vals()
    physical_values = update_physical_values(0, physical_values)
    # compute_numerical_conv_fac(physical_values; return_conv_facs=false, plot_conv_facs=false, print_conv_facs=true)

    conv_facs_analytic = nothing
    param_analytic = nothing

    # # Plot with respect to a_i
    # # var1s = 0:0.1:1
    # # var1_plot_name = L"$a^I$"
    # # xticks = var1s
    # # var1_name = "a_i"

    # # Plot with respecto to a_i and other parameters
    # var1s = [0.1, 0.3, 0.5] # a_i
    # # var2s = [10, 100, 1000, 10000, 100000, 1000000] #t_max
    # # var2s = [15, 50, 100, 200, 300, 500, 700, 1500] # h_atm
    # # var2s = [5, 10, 50, 100, 125, 150, 175, 200] # h_oce
    # # var2s = [20, 200, 2000, 20000] # n_atm
    # # var2s = [5, 50, 500, 5000, 50000] # n_oce
    # # var2s = [10, 50, 100, 125, 150, 175, 200] # L_AI
    # # ai_to_plot = 3
    # # colors = [:blue, :red, :green]
    # # var2s = [0.0007, 0.0008, 0.0009, 0.001, 0.0011, 0.0012, 0.0013] # C_AO
    # # var2s = [0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5] # u_atm and u_oce
    # # var2s = [260, 262, 264, 266, 268, 270, 273, 276] # T_atm_ini
    # # var2s = [270, 271, 272, 273, 274, 275, 276] # T_oce_ini
    # # var2s = [260, 262, 264, 266, 268, 270, 273] # T_ice_ini
    # var2s = 10 .^ LinRange(log10(1), log10(1000), 10)# n_t_atm
    # # var2s = 10 .^ LinRange(log10(1), log10(1000), 10)# n_t_oce

    # xticks = [1, 10, 100, 1000] #var2s
    # var1_name = "a_i"
    # var2_name = "n_t_atm"
    # var2_plot_name = L"$\Delta t^A$"

    # xscale = :log10
    # legend = :topright

    # Run again. See what happens with this better definition of divergence.
    # Evaluating wrt delta z and delta t. Note that delta z should be first argument.
    a_i = 0.1
    delta_z = 10 .^ LinRange(log10(0.1), log10(10), 10)
    n_atms = Int.(round.((physical_values[:h_atm] - physical_values[:z_0numA]) ./ reverse(delta_z)))
    var1s = unique(n_atms)
    physical_values[:Δt_min] = 100
    delta_t = 10 .^ LinRange(log10(1), log10(100), 10)
    n_t_atms = Int.(round.(physical_values[:Δt_min] ./ reverse(delta_t)))
    var2s = unique(n_t_atms)

    # delta_z = 10 .^ LinRange(log10(0.1), log10(10), 10)
    # n_oces = Int.(round.((physical_values[:h_oce] - physical_values[:z_0numO]) ./ delta_z))
    # var1s = unique(n_oces)
    # physical_values[:Δt_min] = 100
    # delta_t = 10 .^ LinRange(log10(1), log10(100), 10)
    # n_t_oces = Int.(round.(physical_values[:Δt_min] ./ delta_t))
    # var2s = unique(n_t_oces)

    yticks = [0.1, 1, 10]#var2s
    xticks = [1, 10, 100]
    var1_name = "n_atm"
    var2_name = "n_t_atm"
    var1_plot_name = L"$\Delta z^A$"
    var2_plot_name = L"$\Delta t^A$"
    xscale = :log10
    yscale = :log10

    # # Computing the convergence factor with slightly different commands.
    # conv_facs_atm, conv_facs_oce, param_analytic, conv_facs_analytic = get_conv_facs_one_variable(physical_values, var1s, var1_name, analytic=true, log_scale=false)
    # conv_facs_atm, conv_facs_oce, param_analytic, conv_facs_analytic = get_conv_facs_two_variables(physical_values, var1s, var2s, var1_name, var2_name, analytic=true, log_scale=(xscale == :log10))
    # # conv_facs_atm, conv_facs_oce, param_analytic, conv_facs_analytic, C_AIs, C_AIs_analytic = get_conv_facs_two_variables(physical_values, var1s, var2s, var1_name, var2_name, analytic=true, log_scale=(xscale == :log10))
    # # conv_facs_atm, conv_facs_oce = get_conv_facs_two_variables(physical_values, var1s, var2s, var1_name, var2_name, analytic=false, log_scale=(xscale == :log10))
    conv_facs_atm, conv_facs_oce = get_conv_facs_two_variables(physical_values, var1s, var2s, var1_name, var2_name, analytic=false, a_i=a_i) # for delta z, delta t dependence.
    # conv_facs_atm, conv_facs_oce = one_iter_update(physical_values, var1s, var2s, var1_name, var2_name)

    conv_facs_numeric = conv_facs_atm

    # Special treatment of n_atm and n_oce.
    if var1_name == "n_atm"
        var1s = (physical_values[:h_atm] - physical_values[:z_0numA]) ./ reverse(var1s)
        conv_facs_numeric = reverse(conv_facs_numeric, dims=1)
    end
    if var2_name == "n_atm"
        var2s = (physical_values[:h_atm] - physical_values[:z_0numA]) ./ reverse(var2s)
        conv_facs_numeric = reverse(conv_facs_numeric, dims=2)
        if !isnothing(conv_facs_analytic) && !isnothing(param_analytic)
            param_analytic = reverse((physical_values[:h_atm] - physical_values[:z_0numA]) ./ param_analytic)
            conv_facs_analytic = reverse(conv_facs_analytic, dims=2)
        end
        xticks = reverse((physical_values[:h_atm] - physical_values[:z_0numA]) ./ xticks)
    end
    if var1_name == "n_oce"
        var1s = (physical_values[:h_oce] - physical_values[:z_0numO]) ./ reverse(var1s)
        conv_facs_numeric = reverse(conv_facs_numeric, dims=1)
    end
    if var2_name == "n_oce"
        var2s = (physical_values[:h_oce] - physical_values[:z_0numO]) ./ reverse(var2s)
        conv_facs_numeric = reverse(conv_facs_numeric, dims=2)
        if !isnothing(conv_facs_analytic) && !isnothing(param_analytic)
            param_analytic = reverse((physical_values[:h_oce] - physical_values[:z_0numO]) ./ param_analytic)
            conv_facs_analytic = reverse(conv_facs_analytic, dims=2)
        end
        xticks = reverse((physical_values[:h_oce] - physical_values[:z_0numO]) ./ xticks)
    end

    # Special treatment for n_t_atm and n_t_oce. Only have as second arguments, so only have to handle this.
    if var2_name == "n_t_atm" || var2_name == "n_t_oce"
        var2s = physical_values[:Δt_min] ./ reverse(var2s)
        conv_facs_numeric = reverse(conv_facs_numeric, dims=2)
        if var1_name == "a_i"
            xticks = physical_values[:Δt_min] ./ reverse(xticks)
        end
        if !isnothing(conv_facs_analytic) && !isnothing(param_analytic)
            param_analytic = reverse(physical_values[:Δt_min] ./ param_analytic)
            conv_facs_analytic = reverse(conv_facs_analytic, dims=2)
        end
    end

    # Special treatment for L_AI
    if var2_name == "L_AI"
        var2s = sort(C_AIs[ai_to_plot, :])
        param_analytic = sort(C_AIs_analytic[ai_to_plot, :])

        rounded_x = round.(var2s, digits=3)
        rounded_and_unique_x = unique(rounded_x)
        # To store the indices of the first occurrences
        first_indices = []

        # Set to track already seen values
        seen_values = Set()

        # Iterate over the vector
        for (i, value) in enumerate(rounded_x)
            if !(value in seen_values)
                push!(first_indices, i)
                push!(seen_values, value)
            end
        end
        x_ticks = var2s[first_indices]
        xticks = (x_ticks, rounded_and_unique_x)
        sorted_C_AI_indices = sortperm(C_AIs[ai_to_plot, :])
        conv_facs_numeric = conv_facs_numeric[ai_to_plot, :][sorted_C_AI_indices]
        var1s = var1s[ai_to_plot]

        sorted_C_AI_indices_analytic = sortperm(C_AIs_analytic[ai_to_plot, :])
        conv_facs_analytic = conv_facs_analytic[ai_to_plot, :][sorted_C_AI_indices_analytic]
    end

    # # Plot wrt ice
    # # plot_wrt_one_param(conv_facs_numeric, var1s, var1_plot_name, conv_facs_analytic=conv_facs_analytic, param_analytic=param_analytic, xticks=xticks, log_conv_fac=true) # For plotting wrt a_i

    # # Plot wrt ice and one other param
    # plot_wrt_a_i_and_one_param(conv_facs_numeric, var1s, var2s, var2_plot_name, conv_facs_analytic=conv_facs_analytic, param_analytic=param_analytic, xscale=xscale, xticks=xticks, scale_text_position=3, text_start_at_bottom=1, legend=legend) # For plotting wrt many params

    # # Plot wrt one parameter (use for C_AI)
    # # plot_wrt_one_param(conv_facs_numeric, var2s, var2_plot_name, conv_facs_analytic=conv_facs_analytic, param_analytic=param_analytic, xticks=xticks, log_conv_fac=false, ylim=:auto, color=colors[ai_to_plot]) # For plotting wrt C_AI

    # # Plot only the numerical convergence factor
    # # plot_wrt_a_i_and_one_param(conv_facs_numeric, var1s, var2s, var2_plot_name, xscale=xscale, xticks=xticks) # For plotting only the numerical. (for instance with regards to delta t)

    # Plot wrt delta t and delta z
    plot_delta_z_delta_t(conv_facs_numeric, var1s, var2s, var1_plot_name, var2_plot_name, xscale=xscale, yscale=yscale, xticks=xticks, yticks=yticks, color=:blue, a_i=a_i)
    # # # Set appropriate backend for plotting
    # plotly()
    # gr()



    # scatter(x, conv_facs_oce[i, :], xticks=x, xscale=:log10, yformatter=:scientific, label="aᴵ=$ai", markersize=5, xlabel=x_axis_label, ylabel="ρᴼ", legend=:topright, color=colors[i])
    # plot!(x, conv_facs_oce[i, :], label="", linewidth=2, color=colors[i])

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
end;

function plot_data(conv_facs_numeric, param_numeric, param_name; conv_facs_analytic=nothing, param_analytic=nothing, plot_numeric=true, xscale=:identity, log_conv_fac=false, label="", color=:black, linestyle=:solid, xticks=:auto, ylim=:auto, legend=:right)
    if !isnothing(conv_facs_analytic) && !isnothing(param_analytic) && plot_numeric
        if log_conv_fac
            y_label = L"$log(̂\hat{\rho})$, $log(\rho)$"
        else
            y_label = L"$\hat{\rho}$, $\rho$"
        end
        # Plot analytical data and a legend
        plot!(param_analytic, conv_facs_analytic, label="analytical" * label, color=color, xscale=xscale, linestyle=linestyle, linewidth=2, legend=legend, xticks=xticks, ylim=ylim)

        # Plot numeric data and a legend
        plot!(param_numeric, conv_facs_numeric, label="numerical" * label, color=color, xscale=xscale, linestyle=linestyle, markershape=:circle, linewidth=2, legend=legend, ylim=ylim)

    elseif !isnothing(conv_facs_analytic) && !isnothing(param_analytic)
        if log_conv_fac
            y_label = L"$log(̂\hat{\rho})$"
        else
            y_label = L"$\hat{\rho}$"
        end
        # Plot analytic data and a legend
        plot!(param_analytic, conv_facs_analytic, label=label, color=color, xscale=xscale, linestyle=linestyle, linewidth=2, legend=legend, xticks=xticks, ylim=ylim)

    elseif plot_numeric
        if log_conv_fac
            y_label = L"$log(̂\rho)$"
        else
            y_label = L"$\rho$"
        end
        # Plot numeric data and a legend
        plot!(param_numeric, conv_facs_numeric, label=label, color=color, xscale=xscale, linestyle=linestyle, markershape=:circle, linewidth=2, legend=legend, xticks=xticks, ylim=ylim)

    end
    xlabel!(param_name)
    ylabel!(y_label)
end

function plot_wrt_a_i_and_one_param(conv_facs_numeric, a_is, param_numeric, param_name; conv_facs_analytic=nothing, param_analytic=nothing, xscale=:identity, xticks=xticks, colors=[:blue, :red, :green], linestyles=[:solid, :dash, :dot], scale_text_position=5, text_start_at_bottom=1, legend=:right)
    """
    Plot conv factor for different a_i and wrt a parameter. Some special treatment
    for divergence has been added if the variable is delta t or delta z
    """
    plot()
    # If there is divergence, to handle the text
    finite_vals = conv_facs_numeric[isfinite.(conv_facs_numeric)]
    if !isnothing(conv_facs_analytic) && !isnothing(param_analytic)
        lowerbound = minimum(conv_facs_analytic)
        upperbound = maximum(conv_facs_analytic)
    elseif !isempty(finite_vals)
        lowerbound = minimum(finite_vals)
        upperbound = maximum(finite_vals)
    else
        lowerbound = 0
        upperbound = 1
    end
    for (i, ai) in enumerate(a_is)
        conv_facs_numeric_i = conv_facs_numeric[i, :]
        param_numeric_not_div = [param_numeric[j] for j in 1:length(conv_facs_numeric_i) if !isinf(conv_facs_numeric_i[j])]
        conv_facs_numeric_i_not_div = [conv_facs_numeric_i[j] for j in 1:length(conv_facs_numeric_i) if !isinf(conv_facs_numeric_i[j])]
        param_numeric_text = [param_numeric[j] for j in 1:length(conv_facs_numeric_i) if isinf(conv_facs_numeric_i[j])]

        # TODO: handle this nicer, dont know how.
        conv_facs_numeric_i_text = [lowerbound + ((i - text_start_at_bottom) * (upperbound - lowerbound) / scale_text_position) for _ in param_numeric_text]
        conv_facs_analytic_i = (isnothing(conv_facs_analytic)) ? nothing : conv_facs_analytic[i, :]
        plot_data(conv_facs_numeric_i_not_div, param_numeric_not_div, param_name, conv_facs_analytic=conv_facs_analytic_i, param_analytic=param_analytic, xscale=xscale, label=", aᴵ = " * "$ai", color=colors[i], linestyle=linestyles[i], xticks=xticks, legend=legend)
        for (k, txt) in enumerate(conv_facs_numeric_i_text)
            annotate!(param_numeric_text[k], txt, text("Divergence", 12, colors[i], :left, rotation=90))
        end
    end
    display(current())
end

function plot_wrt_one_param(conv_facs_numeric, param_numeric, param_name; conv_facs_analytic=nothing, param_analytic=nothing, xticks=nothing, xscale=:identity, log_conv_fac=false, ylim=(-20, -7), color=:black)
    # Plot conv factor as a function of a_i
    plot()
    if log_conv_fac
        conv_facs_analytic = log.(conv_facs_analytic)
        conv_facs_numeric = log.(conv_facs_numeric)
    end
    plot_data(conv_facs_numeric, param_numeric, param_name, conv_facs_analytic=conv_facs_analytic, param_analytic=param_analytic, log_conv_fac=log_conv_fac, xticks=xticks, xscale=xscale, ylim=ylim, color=color)
    display(current())
end

function plot_delta_z_delta_t(conv_facs_numeric, delta_zs, delta_ts, delta_z_name, delta_t_name; xscale=:identity, yscale=:identity, xticks=:auto, yticks=:auto, a_i=0.5, color=:green)
    # Plot convergence factor as a function of deltaz and delta t
    t_values = []
    pre_index = 0
    println(conv_facs_numeric)
    for (i, delta_z) in enumerate(delta_zs)
        first_inf_index = findfirst(isinf, conv_facs_numeric[i, :])

        if first_inf_index === nothing
            push!(t_values, NaN)
        else
            first_t_inf = delta_ts[first_inf_index]
            push!(t_values, first_t_inf)
            if first_inf_index == pre_index
                t_values[i-1] = NaN
            end
        end
        pre_index = first_inf_index
    end

    not_nan_indices = .!isnan.(t_values)
    println(t_values)
    plot(t_values[not_nan_indices], delta_zs[not_nan_indices], linewidth=2, xscale=xscale, yscale=yscale, xticks=xticks, yticks=yticks, xlabel=delta_t_name, ylabel=delta_z_name, label="Divergent regime, aᴵ=$a_i", color=color, fillrange=1e-1, fillalpha=0.3)
    scatter!(t_values[not_nan_indices], delta_zs[not_nan_indices], markershape=:circle, label="", color=color)
    display(current())
end




