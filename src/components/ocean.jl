import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer

struct HeatEquationOcean{P,Y,D,I} <: Interfacer.OceanModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
Interfacer.name(::HeatEquationOcean) = "HeatEquationOcean"


function heat_oce_rhs!(dT, T, cache, t)
    F_sfc = cache.C_AO * (parent(cache.T_air)[1] - T[end])

    ## set boundary conditions
    C3 = CC.Geometry.WVector
    # note: F_sfc is converted to a Cartesian vector in direction 3 (vertical)
    bcs_top = CC.Operators.SetValue(C3(F_sfc))
    bcs_bottom = CC.Operators.SetValue(C3(Float64(0)))

    ## gradient and divergence operators needed for diffusion in tendency calc.
    ᶠgradᵥ = CC.Operators.GradientC2F()
    ᶜdivᵥ = CC.Operators.DivergenceF2C(bottom=bcs_bottom, top=bcs_top)

    @. dT.oce =
        ᶜdivᵥ(cache.k_oce * ᶠgradᵥ(T.oce)) /
        (cache.ρ_oce * cache.c_oce)
end

function ocean_init(stepping, ics, space, cache)
    Δt = Float64(stepping.Δt_min) / stepping.nsteps_oce
    saveat = Float64(stepping.Δt_coupler)

    ode_function = CTS.ClimaODEFunction((T_exp!)=heat_oce_rhs!)

    problem = SciMLBase.ODEProblem(ode_function, ics, stepping.timerange, cache)
    integrator = SciMLBase.init(
        problem,
        stepping.odesolver,
        dt=Δt,
        saveat=saveat,
        adaptive=false,
    )

    sim = HeatEquationOcean(cache, ics, space, integrator)
    return sim
end

get_field(sim::HeatEquationOcean, ::Val{:T_oce_sfc}) = sim.integrator.u[end]
function update_field!(sim::HeatEquationOcean, ::Val{:T_atm_sfc}, field)
    parent(sim.integrator.p.T_air) .= field
end


Interfacer.step!(sim::HeatEquationOcean, t) =
    Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.reinit!(sim::HeatEquationOcean) = Interfacer.reinit!(sim.integrator)

Checkpointer.get_model_prog_state(sim::HeatEquationOcean) = sim.integrator.u
