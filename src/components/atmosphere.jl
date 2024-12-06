import ClimaCore as CC
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer

struct HeatEquationAtmos{P, Y, D, I} <: Interfacer.AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
Interfacer.name(::HeatEquationAtmos) = "HeatEquationAtmos"

function heat_atm_rhs!(dT, T, cache, t)
    F_sfc = cache.params.C_AO * (parent(T)[1] - cache.T_sfc)

    ## set boundary conditions
    C3 = CC.Geometry.WVector
    # note: F_sfc is converted to a Cartesian vector in direction 3 (vertical)
    bcs_bottom = CC.Operators.SetValue(C3(F_sfc))
    bcs_top = CC.Operators.SetValue(C3(Float64(0)))

    ## gradient and divergence operators needed for diffusion in tendency calc.
    ᶠgradᵥ = CC.Operators.GradientC2F()
    ᶜdivᵥ = CC.Operators.DivergenceF2C(bottom = bcs_bottom, top = bcs_top)

    @. dT.atm = ᶜdivᵥ(cache.params.k_atm * ᶠgradᵥ(T.atm)) / (cache.params.ρ_atm * cache.params.c_atm)
end

function atmos_init(
    stepping,
    ics,
    space,
    parameters,
)
    Δt = Float64(stepping.Δt_min)

    cache = (
        params = parameters,
        F_turb_energy = CC.Fields.zeros(space),
        F_radiative = CC.Fields.zeros(space),
        q_sfc = CC.Fields.zeros(space),
        ρ_sfc = CC.Fields.zeros(space),
        T_sfc = Float64(0),
    )

    ode_function = CTS.ClimaODEFunction(T_exp! = heat_atm_rhs!)
    problem = SciMLBase.ODEProblem(ode_function, ics, stepping.timerange, cache)
    integrator = SciMLBase.init(problem, stepping.odesolver, dt = Δt, saveat = Float64(stepping.Δt_coupler), adaptive = false)

    sim = HeatEquationAtmos(parameters, ics, space, integrator)
    return sim
end

Checkpointer.get_model_prog_state(sim::HeatEquationAtmos) = sim.integrator.u

Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:radiative_energy_flux_toa}) = nothing
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:energy}) = nothing
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:air_density}) = sim.integrator.p.ρ_atm
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:air_temperature}) = sim.integrator.u[1]
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:liquid_precipitation}) = sim.integrator.u
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:snow_precipitation}) = sim.integrator.u
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:radiative_energy_flux_sfc}) = sim.integrator.u
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:turbulent_energy_flux}) = nothing
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:turbulent_moisture_flux}) = nothing
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:thermo_state_int}) = nothing
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:water}) = nothing
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:height_int}) = nothing
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:height_sfc}) = nothing
Interfacer.get_field(sim::HeatEquationAtmos, ::Val{:uv_int}) = nothing


Interfacer.update_field!(sim::HeatEquationAtmos, ::Val{:turbulent_fluxes}, fields) = nothing
Interfacer.update_field!(sim::HeatEquationAtmos, ::Val{:surface_direct_albedo}, fields) = nothing

# extensions required by FieldExchanger
Interfacer.step!(sim::HeatEquationAtmos, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.reinit!(sim::HeatEquationAtmos) = Interfacer.reinit!(sim.integrator)


FluxCalculator.calculate_surface_air_density(sim::HeatEquationAtmos, T_S::CC.Fields.Field) = T_S

function FluxCalculator.atmos_turbulent_fluxes_most!(atmos_sim::HeatEquationAtmos, csf)
    p = atmos_sim.integrator.p
    p.turbulent_energy_flux = p.transfer_coefficient * (csf.T_atm - csf.T_oce)
end

