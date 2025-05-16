import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer

export HeatEquationAtmos, heat_atm_rhs!, atmos_init, get_field, update_field!

struct HeatEquationAtmos{P,Y,D,I} <: Interfacer.AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

Interfacer.name(::HeatEquationAtmos) = "HeatEquationAtmos"

function heat_atm_rhs!(dT, T, cache, t)
    if cache.boundary_mapping == "mean"
        F_sfc = (
            cache.a_I *
            cache.C_AI *
            (T[1] - parent(cache.T_Is)[1]) +
            (1 - cache.a_I) *
            cache.C_AO *
            (T[1] - parent(cache.T_O)[1])
        )
    else
        index = argmin(abs.(parent(CC.Fields.coordinate_field(cache.T_O)) .- t))
        F_sfc = (
            cache.a_I *
            cache.C_AI *
            (T[1] - parent(cache.T_Is)[index]) +
            (1 - cache.a_I) *
            cache.C_AO *
            (T[1] - parent(cache.T_O)[index])
        )
    end
    # set boundary conditions
    C3 = CC.Geometry.WVector
    # note: F_sfc is converted to a Cartesian vector in direction 3 (vertical)
    bcs_bottom = CC.Operators.SetValue(C3(F_sfc))
    bcs_top = CC.Operators.SetValue(C3(Float64(0)))

    ## gradient and divergence operators needed for diffusion in tendency calc.
    ᶠgradᵥ = CC.Operators.GradientC2F()
    ᶜdivᵥ = CC.Operators.DivergenceF2C(bottom=bcs_bottom, top=bcs_top)

    @. dT.data = ᶜdivᵥ(cache.k_A * ᶠgradᵥ(T.data)) / (cache.ρ_A * cache.c_A)
end

function atmos_init(stepping, ics, space, cache)
    Δt = Float64(stepping.Δt_min) / stepping.nsteps_atm
    saveat = stepping.timerange[1]:stepping.Δt_min:stepping.timerange[end]

    ode_function = CTS.ClimaODEFunction((T_exp!)=heat_atm_rhs!)
    problem = SciMLBase.ODEProblem(ode_function, ics, stepping.timerange, cache)
    integrator = SciMLBase.init(
        problem,
        stepping.odesolver,
        dt=Δt,
        saveat=saveat,
        adaptive=false,
    )

    sim = HeatEquationAtmos(cache, ics, space, integrator)
    return sim
end

Checkpointer.get_model_prog_state(sim::HeatEquationAtmos) = sim.integrator.u

function Interfacer.step!(sim::HeatEquationAtmos, t)
    Interfacer.step!(sim.integrator, t - sim.integrator.t)
    check_stability(sim.integrator.u, sim.params.stable_range)
end

Interfacer.reinit!(sim::HeatEquationAtmos) = Interfacer.reinit!(sim.integrator)

function get_field(sim::HeatEquationAtmos, ::Val{:T_atm_sfc})
    return vec([fieldvec[1] for fieldvec in sim.integrator.sol.u])
end

function update_field!(sim::HeatEquationAtmos, T_O, T_Is)
    if sim.params.boundary_mapping == "mean"
        T_O = vec([mean(T_O)])
        T_Is = vec([mean(T_Is)])
    end
    parent(sim.integrator.p.T_O) .= T_O
    parent(sim.integrator.p.T_Is) .= T_Is
end
