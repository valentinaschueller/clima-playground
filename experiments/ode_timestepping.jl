import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Interfacer
import LinearAlgebra
import SciMLBase
import Plots
import ClimaCore.MatrixFields: @name, FieldMatrixWithSolver, FieldMatrix, DiagonalMatrixRow

struct Land{P,Y,D,I} <: Interfacer.LandModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

Interfacer.step!(sim::Land, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t)

function ∑sfc_flx!(dT, T, p, t)
    @. dT = p.λ
end

function Wfact(W, Y, p, dtγ, t)
    @. W.matrix[@name(data), @name(data)] = dtγ * p.λ * (LinearAlgebra.I,) - (LinearAlgebra.I,)
end

T_0 = 1.0
t_0 = 0.0
t_end = 1.0
Δt = 0.1
p = (λ=-0.5,)

context = CC.ClimaComms.context()
point_space = CC.Spaces.PointSpace(context, CC.Geometry.ZPoint(0.0))
ics = CC.Fields.FieldVector(data=CC.Fields.ones(point_space) .* T_0)

timestepping = :implicit
if timestepping == :implicit
    jacobian = CC.MatrixFields.FieldMatrix(
        (@name(data), @name(data)) => similar(ics.data, CC.MatrixFields.DiagonalMatrixRow{Float64}),
    )
    T_imp! = SciMLBase.ODEFunction(∑sfc_flx!; jac_prototype=FieldMatrixWithSolver(jacobian, ics), Wfact=Wfact)
    odesolver = CTS.IMEXAlgorithm(CTS.ARS111(), CTS.NewtonsMethod())
    ode_function = CTS.ClimaODEFunction((T_imp!)=T_imp!)
else
    odesolver = CTS.ExplicitAlgorithm(CTS.RK4())
    ode_function = CTS.ClimaODEFunction((T_exp!)=∑sfc_flx!)
end

problem = SciMLBase.ODEProblem(ode_function, ics, (t_0, t_end), p)

integrator = SciMLBase.init(
    problem,
    odesolver,
    dt=Δt,
    saveat=t_0:Δt:t_end,
    adaptive=false,
)
sim = Land(p, ics, point_space, integrator)
Interfacer.step!(sim, t_end)

SciMLBase.solve!(integrator);

solution = [parent(fieldvec.data)[1] for fieldvec in integrator.sol.u]
Plots.plot(integrator.sol.t, solution)
