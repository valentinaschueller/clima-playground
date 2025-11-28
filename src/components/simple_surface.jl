import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, Interfacer, Utilities, FluxCalculator
import ClimaDiagnostics as CD
import ClimaCore.MatrixFields: @name, FieldMatrixWithSolver, FieldMatrix
import ClimaAtmos as CA

export SimpleSeaIce, simple_surface_init

struct SimpleSeaIce{P,Y,D,I} <: Interfacer.SeaIceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

function sfc_rhs!(dh, h, p, t)
    dh.data = 0.0
    return
end


function simple_surface_init(p::AbstractDict, output_dir)
    FT = p["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
    context = Utilities.get_comms_context(Dict("device" => "auto"))
    space = CC.Spaces.PointSpace(context, CC.Geometry.ZPoint(FT(0.0)))
    field_h_I = CC.Fields.ones(space) .* p["h_I_ini"]
    ics = CC.Fields.FieldVector(data=field_h_I)

    t_end = CA.time_to_seconds(p["t_end"])

    ode_function = CTS.ClimaODEFunction((T_exp!)=sfc_rhs!)
    problem = SciMLBase.ODEProblem(ode_function, ics, (FT(0), t_end), p)
    Δt = CA.time_to_seconds(p["dt"])

    saveat = 0.0:Δt:t_end

    integrator = SciMLBase.init(
        problem,
        CTS.ExplicitAlgorithm(CTS.RK4()),
        dt=Δt,
        saveat=saveat,
        adaptive=false,
    )
    sim = SimpleSeaIce(p, ics, space, integrator)
    return sim
end

Checkpointer.get_model_prog_state(sim::SimpleSeaIce) = sim.integrator.u

function Interfacer.step!(sim::SimpleSeaIce, t)
    Interfacer.step!(sim.integrator, t - sim.integrator.t)
end

function Interfacer.get_field(sim::SimpleSeaIce, ::Val{:h_I})
    return sim.integrator.p.h_I
end

function Interfacer.get_field(sim::SimpleSeaIce, ::Val{:area_fraction})
    return 1.0
end

function Interfacer.get_field(sim::SimpleSeaIce, ::Val{:ice_concentration})
    return 1.0
end

function Interfacer.get_field(sim::SimpleSeaIce, ::Val{:emissivity})
    return sim.integrator.p["eps_I"]
end

function Interfacer.get_field(sim::SimpleSeaIce, ::Val{:surface_temperature})
    return sim.integrator.p["T_sfc"]
end

function Interfacer.get_field(sim::SimpleSeaIce, ::Val{:surface_diffuse_albedo})
    return sim.integrator.p["alb_I"]
end

function Interfacer.get_field(sim::SimpleSeaIce, ::Val{:surface_direct_albedo})
    return sim.integrator.p["alb_I"]
end


Interfacer.update_field!(sim::SimpleSeaIce, ::Union{
        Val{:area_fraction},
        Val{:SW_d},
        Val{:LW_d},
        Val{:snow_precipitation},
        Val{:liquid_precipitation},
        Val{:turbulent_energy_flux},
        Val{:turbulent_moisture_flux}
    }, field) = nothing

Interfacer.get_field(sim::SimpleSeaIce, ::Val{:roughness_momentum}) = sim.integrator.p["z_AO"]
Interfacer.get_field(sim::SimpleSeaIce, ::Val{:roughness_buoyancy}) = sim.integrator.p["z_AO"]

FluxCalculator.update_turbulent_fluxes!(sim::SimpleSeaIce, fields::NamedTuple) = nothing

function Interfacer.add_coupler_fields!(coupler_field_names, ::SimpleSeaIce)
    coupler_fields = [
        :surface_temperature, 
        :surface_diffuse_albedo, 
        :surface_direct_albedo, 
        :emissivity, 
        :area_fraction, 
        :ice_concentration,
    ]
    push!(coupler_field_names, coupler_fields...)
end


