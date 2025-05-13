import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer

export HeatEquationOcean, heat_oce_rhs!, ocean_init, get_field, update_field!

struct HeatEquationOcean{P,Y,D,I} <: Interfacer.OceanModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
Interfacer.name(::HeatEquationOcean) = "HeatEquationOcean"

function heat_oce_rhs!(dT, T, cache, t)
    if cache.boundary_mapping == "mean"
        F_sfc = (
            cache.a_I *
            cache.C_IO *
            (cache.T_Ib - T[end]) +
            (1 - cache.a_I) *
            cache.C_AO *
            (parent(cache.T_A)[1] - T[end])
        )
    else
        index = argmin(abs.(parent(CC.Fields.coordinate_field(cache.T_A)) .- t))
        F_sfc = (
            cache.a_I *
            cache.C_IO *
            (cache.T_Ib - T[end]) +
            (1 - cache.a_I) *
            cache.C_AO *
            (parent(cache.T_A)[index] - T[end])
        )
    end

    ## set boundary conditions
    C3 = CC.Geometry.WVector
    # note: F_sfc is converted to a Cartesian vector in direction 3 (vertical)
    bcs_top = CC.Operators.SetValue(C3(F_sfc))
    bcs_bottom = CC.Operators.SetValue(C3(Float64(0)))

    ## gradient and divergence operators needed for diffusion in tendency calc.
    ᶠgradᵥ = CC.Operators.GradientC2F()
    ᶜdivᵥ = CC.Operators.DivergenceF2C(bottom=bcs_bottom, top=bcs_top)

    @. dT.data = ᶜdivᵥ(cache.k_O * ᶠgradᵥ(T.data)) / (cache.ρ_O * cache.c_O)
end

function ocean_init(stepping, ics, space, cache)
    Δt = Float64(stepping.Δt_min) / stepping.nsteps_oce
    saveat = stepping.timerange[1]:stepping.Δt_min:stepping.timerange[end]

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

Checkpointer.get_model_prog_state(sim::HeatEquationOcean) = sim.integrator.u

function Interfacer.step!(sim::HeatEquationOcean, t)
    Interfacer.step!(sim.integrator, t - sim.integrator.t)
    if !is_stable(sim.integrator.u, sim.params.stable_range)
        throw(UnstableError())
    end
end


Interfacer.reinit!(sim::HeatEquationOcean) = Interfacer.reinit!(sim.integrator)

get_field(sim::HeatEquationOcean, ::Val{:T_oce_sfc}) = sim.integrator.u[end]

function update_field!(sim::HeatEquationOcean, field_1)
    if sim.params.boundary_mapping == "mean"
        parent(sim.integrator.p.T_A)[1] = field_1
    else
        parent(sim.integrator.p.T_A) .= field_1
    end
end

