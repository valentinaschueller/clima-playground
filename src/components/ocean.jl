import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, FluxCalculator, Interfacer, Utilities

struct HeatEquationOcean{P,Y,D,I} <: Interfacer.OceanModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
Interfacer.name(::HeatEquationOcean) = "HeatEquationOcean"


function heat_oce_rhs!(dT, T, cache, t)
    F_sfc = cache.params.C_AO * (cache.T_air - last(parent(T)))

    ## set boundary conditions
    C3 = CC.Geometry.WVector
    # note: F_sfc is converted to a Cartesian vector in direction 3 (vertical)
    bcs_top = CC.Operators.SetValue(C3(F_sfc))
    bcs_bottom = CC.Operators.SetValue(C3(Float64(0)))

    ## gradient and divergence operators needed for diffusion in tendency calc.
    ᶠgradᵥ = CC.Operators.GradientC2F()
    ᶜdivᵥ = CC.Operators.DivergenceF2C(bottom=bcs_bottom, top=bcs_top)

    @. dT.oce =
        ᶜdivᵥ(cache.params.k_oce * ᶠgradᵥ(T.oce)) /
        (cache.params.ρ_oce * cache.params.c_oce)
end

function ocean_init(stepping, ics, space, parameters)
    Δt = Float64(stepping.Δt_min)
    saveat = Float64(stepping.Δt_coupler)

    cache = (params=parameters, area_fraction=1 - parameters.a_i, T_air=Float64(0))

    ode_function = CTS.ClimaODEFunction((T_exp!)=heat_oce_rhs!)

    problem = SciMLBase.ODEProblem(ode_function, ics, stepping.timerange, cache)
    integrator = SciMLBase.init(
        problem,
        stepping.odesolver,
        dt=Δt,
        saveat=saveat,
        adaptive=false,
    )

    sim = HeatEquationOcean(parameters, ics, space, integrator)
    return sim
end

Interfacer.get_field(sim::HeatEquationOcean, ::Val{:air_density}) = nothing
Interfacer.get_field(sim::HeatEquationOcean, ::Val{:area_fraction}) =
    sim.integrator.p.area_fraction
Interfacer.get_field(sim::HeatEquationOcean, ::Val{:beta}) = nothing
Interfacer.get_field(sim::HeatEquationOcean, ::Val{:roughness_buoyancy}) = nothing
Interfacer.get_field(sim::HeatEquationOcean, ::Val{:roughness_momentum}) = nothing
Interfacer.get_field(sim::HeatEquationOcean, ::Val{:surface_direct_albedo}) = Float64(0)
Interfacer.get_field(sim::HeatEquationOcean, ::Val{:surface_diffuse_albedo}) = Float64(0)
Interfacer.get_field(sim::HeatEquationOcean, ::Val{:surface_humidity}) = nothing
Interfacer.get_field(sim::HeatEquationOcean, ::Val{:surface_temperature}) =
    sim.integrator.u[end]
Interfacer.get_field(sim::HeatEquationOcean, ::Val{:water}) = nothing


Interfacer.get_field(sim::HeatEquationOcean, ::Val{:energy}) =
    sim.integrator.p.params.ρ .* sim.integrator.p.params.c .* sim.integrator.u.oce .*
    sim.integrator.p.params.h

Interfacer.update_field!(
    sim::HeatEquationOcean,
    ::Val{:area_fraction},
    field::CC.Fields.Field,
) = nothing
Interfacer.update_field!(
    sim::HeatEquationOcean,
    ::Val{:air_density},
    field::CC.Fields.Field,
) = nothing
Interfacer.update_field!(
    sim::HeatEquationOcean,
    ::Val{:radiative_energy_flux_sfc},
    field::CC.Fields.Field,
) = nothing
Interfacer.update_field!(
    sim::HeatEquationOcean,
    ::Val{:turbulent_energy_flux},
    field::CC.Fields.Field,
) = Float64(0)
Interfacer.update_field!(
    sim::HeatEquationOcean,
    ::Val{:surface_direct_albedo},
    field::CC.Fields.Field,
) = Float64(0)
Interfacer.update_field!(
    sim::HeatEquationOcean,
    ::Val{:surface_diffuse_albedo},
    field::CC.Fields.Field,
) = nothing

Interfacer.step!(sim::HeatEquationOcean, t) =
    Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.reinit!(sim::HeatEquationOcean) = Interfacer.reinit!(sim.integrator)

# extensions required by FluxCalculator (partitioned fluxes)
FluxCalculator.update_turbulent_fluxes!(sim::HeatEquationOcean, fields::NamedTuple) =
    nothing

function Checkpointer.get_model_prog_state(sim::HeatEquationOcean)
    return sim.integrator.u
end
