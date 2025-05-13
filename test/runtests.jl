import Test: @test
import Statistics
using clima_playground

cs, _ = coupled_heat_equations(iterations=10, Δt_min=10, t_max=1000, Δt_cpl=1000)
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.sol.u[end]), 267.0199795729795)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.9999761867679)

cs, _ = coupled_heat_equations(iterations=10, parallel=true, Δt_min=10, t_max=1000, Δt_cpl=1000)
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.sol.u[end]), 267.0199795729795)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.9999761867679)

cs, _ = coupled_heat_equations(iterations=10, Δt_min=10, t_max=1000, Δt_cpl=1000, boundary_mapping="cit")
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.sol.u[end]), 267.0199795729795)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.9999761867679)

cs, _ = coupled_heat_equations(iterations=10, parallel=true, Δt_min=10, t_max=1000, Δt_cpl=1000, boundary_mapping="cit")
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.sol.u[end]), 267.0199795729795)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.9999761867679)

cs, _ = coupled_heat_equations(iterations=10, a_I=0.5, Δt_min=10, t_max=1000, Δt_cpl=1000, T_Ib=270.0)
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.sol.u[end]), 267.0157868000668)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.9816379915888)
@test isapprox(Statistics.mean(cs.model_sims.ice_sim.integrator.sol.u[end]), 270.0)

cs, _ = coupled_heat_equations(iterations=10, a_I=1.0, Δt_min=10, T_I_ini=250.0, t_max=1000, Δt_cpl=1000, T_Ib=250.0)
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.sol.u[end]), 266.9138392636568)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.5827944417085)
@test isapprox(Statistics.mean(cs.model_sims.ice_sim.integrator.sol.u[end]), 250.0)

cs, _ = coupled_heat_equations()
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.sol.u[end]), 267.021160864323)
@test isapprox(maximum(cs.model_sims.atmos_sim.integrator.sol.u[end]), 270.98203813664367)
@test isapprox(minimum(cs.model_sims.atmos_sim.integrator.sol.u[end]), 267.0)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.99996464720743)
@test isapprox(maximum(cs.model_sims.ocean_sim.integrator.sol.u[end]), 271.0)
@test isapprox(minimum(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.998233166346)

p = SimulationParameters()
p.C_H_AO = compute_C_H_AO(p)
p.C_H_AI = compute_C_H_AI(p)
restore_physical_values!(p)
@test isapprox(compute_ϱ_analytical(p), 0.10729333488914307)