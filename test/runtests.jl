import Test: @test
import Statistics
using clima_playground

cs, _ = coupled_heat_equations(iterations=10, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000))
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.sol.u[end]), 267.0199795729795)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.9999761867679)

cs, _ = coupled_heat_equations(iterations=10, parallel=true, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000))
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.sol.u[end]), 267.0199795729795)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.9999761867679)

cs, _ = coupled_heat_equations(iterations=10, boundary_mapping="cit", params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000))
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.sol.u[end]), 267.0199795729795)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.9999761867679)

cs, _ = coupled_heat_equations(iterations=10, boundary_mapping="cit", parallel=true, params=Dict(:Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000))
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.sol.u[end]), 267.0199795729795)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.9999761867679)

cs, _ = coupled_heat_equations(iterations=10, params=Dict(:a_i => 0.5, :Δt_min => 10, :t_max => 1000, :Δt_cpl => 1000))
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.sol.u[end]), 267.0157868000668)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.sol.u[end]), 270.9816379915888)
@test isapprox(Statistics.mean(cs.model_sims.ice_sim.integrator.sol.u[end]), 270.0)

cs, _ = coupled_heat_equations(iterations=10, params=Dict(:a_i => 1.0, :Δt_min => 10, :T_ice_ini => 250.0, :t_max => 1000, :Δt_cpl => 1000))
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

@test isapprox(compute_ρ_analytical(define_realistic_vals()), 0.10729333488914307)