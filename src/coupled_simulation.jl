import Dates
import SciMLBase
import ClimaAtmos as CA
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Interfacer, Utilities
import YAML

export get_coupled_sim, get_odesolver


function get_odesolver(::Val{:implicit})
    return CTS.IMEXAlgorithm(CTS.ARS111(), CTS.NewtonsMethod())
end

function get_odesolver(::Val{:explicit})
    return CTS.ExplicitAlgorithm(CTS.RK4())
end

function get_coupled_sim(p::SimulationParameters)
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
    config = CA.AtmosConfig("experiments/config.yaml")
    atmos_sim = ClimaAtmosSimulation(config)
    thermo_params = get_thermo_params(atmos_sim)

    ocean_sim = ocean_init(odesolver, p, output_dir)
    ice_sim = ice_init(odesolver, p, output_dir)

    boundary_space = ice_sim.domain

    start_date = Dates.DateTime("19790301", Dates.dateformat"yyyymmdd")

    model_sims = (atmos_sim=atmos_sim, ocean_sim=ocean_sim, ice_sim=ice_sim)

    coupler_field_names = []
    for sim in model_sims
        Interfacer.add_coupler_fields!(coupler_field_names, sim)
    end
    coupler_fields = Interfacer.init_coupler_fields(Float64, coupler_field_names, boundary_space)

    tspan = (p.t_0, p.t_max)
    cs = Interfacer.CoupledSimulation{Float64}(
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
        thermo_params, # thermo_params
        nothing, # diagnostic_handlers
    )
    return cs
end
