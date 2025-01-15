import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer
using CSV
using DataFrames

struct HeatEquationAtmos{P,Y,D,I} <: Interfacer.AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

Interfacer.name(::HeatEquationAtmos) = "HeatEquationAtmos"

function heat_atm_rhs!(dT, T, cache, t) # Apparently this is called 4 times... per time point
    F_sfc = (cache.a_i * cache.C_AI * (T[1] - parent(cache.T_ice)[1]) + (1 - cache.a_i) * cache.C_AO * (T[1] - parent(cache.T_sfc)[1])) / cache.k_atm # I say we should divide by k^A here?
    # set boundary conditions
    C3 = CC.Geometry.WVector
    # note: F_sfc is converted to a Cartesian vector in direction 3 (vertical)
    bcs_bottom = CC.Operators.SetValue(C3(F_sfc))
    bcs_top = CC.Operators.SetValue(C3(Float64(0)))

    ## gradient and divergence operators needed for diffusion in tendency calc.
    ᶠgradᵥ = CC.Operators.GradientC2F()
    ᶜdivᵥ = CC.Operators.DivergenceF2C(bottom=bcs_bottom, top=bcs_top) # Do we not miss k^A here for the bcs?

    @. dT.atm =
        ᶜdivᵥ(cache.k_atm * ᶠgradᵥ(T.atm)) /
        (cache.ρ_atm * cache.c_atm)
end

function atmos_init(stepping, ics, space, cache)
    Δt = Float64(stepping.Δt_min) / stepping.nsteps_atm

    ode_function = CTS.ClimaODEFunction((T_exp!)=heat_atm_rhs!)
    problem = SciMLBase.ODEProblem(ode_function, ics, stepping.timerange, cache)
    integrator = SciMLBase.init(
        problem,
        stepping.odesolver,
        dt=Δt,
        saveat=Float64(stepping.Δt_min), # Change here?
        adaptive=false,
    )

    sim = HeatEquationAtmos(cache, ics, space, integrator)
    return sim
end

Checkpointer.get_model_prog_state(sim::HeatEquationAtmos) = sim.integrator.u

Interfacer.step!(sim::HeatEquationAtmos, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t)
Interfacer.reinit!(sim::HeatEquationAtmos) = Interfacer.reinit!(sim.integrator)

get_field(sim::HeatEquationAtmos, ::Val{:T_atm_sfc}) = sim.integrator.u[1]
function update_field!(sim::HeatEquationAtmos, ::Val{:T_oce_sfc}, field_1, ::Val{:T_ice}, field_2)
    parent(sim.integrator.p.T_sfc)[1] = field_1
    parent(sim.integrator.p.T_ice)[1] = field_2
end






