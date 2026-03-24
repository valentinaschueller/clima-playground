import Dates
import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaCoupler: FieldExchanger, Interfacer

export get_coupled_sim, get_odesolver

function get_odesolver(::Val{:implicit})
    return CTS.IMEXAlgorithm(CTS.ARS111(), CTS.NewtonsMethod())
end

function get_odesolver(::Val{:explicit})
    return CTS.ExplicitAlgorithm(CTS.RK4())
end

function get_coupled_sim(p::SimulationParameters)
    FT = eltype(p)
    output_dir = "output"
    rm(output_dir, recursive=true, force=true)
    mkpath(output_dir)
    dir_paths = (
        output=output_dir,
        artifacts=output_dir,
        regrid=output_dir,
        checkpoints=output_dir,
    )

    p.T_A = p.T_A_ini
    p.T_O = p.T_O_ini
    p.F_AO = p.C_AO * (p.T_A - p.T_O)

    if p.ice_model_type != :constant
        @info("Determine initial ice surface temperature from SEB.")
        p.T_I_ini = T_Is(p)
    end
    p.T_Is = p.T_I_ini

    if p.ice_model_type == :constant
        initial_values = [p.T_A, p.T_O, p.T_Is, p.T_Ib]
        min_value = minimum(initial_values)
        max_value = maximum(initial_values)
        p.stable_range = (min_value - eps(min_value), max_value + eps(max_value))
    end

    odesolver = get_odesolver(Val(p.timestepping))
    atmos_sim = atmos_init(odesolver, p, output_dir)
    ocean_sim = ocean_init(odesolver, p, output_dir)
    ice_sim = ice_init(odesolver, p, output_dir)

    boundary_space = ice_sim.domain

    start_date = Dates.DateTime("19790301", Dates.dateformat"yyyymmdd")

    model_sims = (atmos_sim=atmos_sim, ocean_sim=ocean_sim, ice_sim=ice_sim)

    coupler_field_names = []
    for sim in model_sims
        Interfacer.add_coupler_fields!(coupler_field_names, sim)
    end
    coupler_fields = Interfacer.init_coupler_fields(FT, coupler_field_names, boundary_space)

    tspan = (p.t_0, p.t_max)
    cs = Interfacer.CoupledSimulation{FT}(
        Ref(start_date),
        coupler_fields,
        nothing, # conservation checks
        tspan,
        p.Δt_cpl,
        Ref(tspan[begin]),
        Ref(-1),
        model_sims,
        (;), # callbacks
        dir_paths,
        nothing, # thermo_params
        nothing, # diagnostic_handlers
    )
    return cs
end

"""Sets u0 = u(t), t0 = t, and empties u."""
function reinit!(cs::Interfacer.CoupledSimulation, t)
    for sim in cs.model_sims
        SciMLBase.reinit!(sim.integrator, sim.integrator.u, t0=t)
    end
end

function update_atmos_values!(cs)
    bound_ocean_vals = Interfacer.get_field(cs.model_sims.ocean_sim, Val(:T_oce_sfc))
    ice_T = Interfacer.get_field(cs.model_sims.ice_sim, Val(:T_ice))
    Interfacer.update_field!(cs.model_sims.atmos_sim, bound_ocean_vals, ice_T)
    return bound_ocean_vals
end

function update_ocean_values!(cs)
    bound_atmos_vals = Interfacer.get_field(cs.model_sims.atmos_sim, Val(:F_AO))
    Interfacer.update_field!(cs.model_sims.ocean_sim, bound_atmos_vals)
    return bound_atmos_vals
end

function update_ice_values!(cs)
    bound_atmos_vals = Interfacer.get_field(cs.model_sims.atmos_sim, Val(:T_atm_sfc))
    bound_ocean_vals = Interfacer.get_field(cs.model_sims.ocean_sim, Val(:T_oce_sfc))
    Interfacer.update_field!(cs.model_sims.ice_sim, bound_atmos_vals, bound_ocean_vals)
    return bound_atmos_vals, bound_ocean_vals
end

function advance_simulation!(cs::Interfacer.CoupledSimulation, t_end::FT, parallel::Bool) where {FT}
    if parallel
        FieldExchanger.step_model_sims!(cs.model_sims, t_end)
        update_atmos_values!(cs)
        update_ocean_values!(cs)
        bound_atmos_vals, bound_ocean_vals = update_ice_values!(cs)
    else
        Interfacer.step!(cs.model_sims.ice_sim, t_end)
        Interfacer.step!(cs.model_sims.ocean_sim, t_end)
        update_atmos_values!(cs)
        Interfacer.step!(cs.model_sims.atmos_sim, t_end)
        update_ocean_values!(cs)
        bound_atmos_vals, bound_ocean_vals = update_ice_values!(cs)
    end
    return bound_atmos_vals, bound_ocean_vals
end


function solve_coupler!(
    cs::Interfacer.CoupledSimulation;
    iterations::Int=1,
    parallel::Bool=false,
)
    (; Δt_cpl, tspan) = cs

    @info("Starting coupling loop")

    ϱ_A = nothing
    ϱ_O = nothing

    t0 = tspan[begin]

    for t = t0:Δt_cpl:(tspan[end]-Δt_cpl)
        @info("Current time: $t")
        cs.t[] = t
        reinit!(cs, t)

        iter = 1
        atmos_vals_list = []
        ocean_vals_list = []
        T_A_Γ = nothing
        T_O_Γ = nothing

        Checkpointer.checkpoint_sims(cs)

        for iter = 1:iterations
            @info("Current iter: $(iter)")
            if iter > 1
                Checkpointer.restart!(cs, cs.dirs.checkpoints, Int(cs.t[]))
                reinit!(cs, t)
            end

            T_A_Γ_old = T_A_Γ
            T_O_Γ_old = T_O_Γ

            try
                T_A_Γ, T_O_Γ = advance_simulation!(cs, t + Δt_cpl, parallel)
            catch err
                if isa(err, UnstableError)
                    @warn("Unstable simulation!")
                    return NaN, NaN
                end
                rethrow()
            end

            push!(atmos_vals_list, T_A_Γ)
            push!(ocean_vals_list, T_O_Γ)

            if has_converged(T_A_Γ, T_A_Γ_old, T_O_Γ, T_O_Γ_old)
                @info "Termination criterion satisfied!"
            end
        end
        if iterations > 1
            ϱ_A = compute_ϱ_numerical(atmos_vals_list, parallel)
            ϱ_O = compute_ϱ_numerical(ocean_vals_list, parallel)
        end
    end
    cs.t[] = tspan[end]
    return ϱ_A, ϱ_O
end
