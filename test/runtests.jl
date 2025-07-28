import Test: @test
import Statistics
using clima_playground

cs, _ = coupled_heat_equations(iterations=10, Δt_min=10, t_max=1000, Δt_cpl=1000, timestepping=:explicit)
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.u), 267.0199795729795)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.u), 270.9999761867679)

cs, _ = coupled_heat_equations(iterations=10, parallel=true, Δt_min=10, t_max=1000, Δt_cpl=1000, timestepping=:explicit)
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.u), 267.0199795729795)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.u), 270.9999761867679)

cs, _ = coupled_heat_equations(iterations=10, a_I=0.5, Δt_min=10, t_max=1000, Δt_cpl=1000, T_Ib=270.0, timestepping=:explicit)
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.u), 267.0157868000668)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.u), 270.9816379915888)
@test isapprox(Statistics.mean(get_field(cs.model_sims.ice_sim, Val(:T_ice))), 270.0)

cs, _ = coupled_heat_equations(iterations=10, a_I=1.0, Δt_min=10, T_I_ini=250.0, t_max=1000, Δt_cpl=1000, T_Ib=250.0, timestepping=:explicit)
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.u), 266.9138392636568)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.u), 270.5827944417085)
T_Is = get_field(cs.model_sims.ice_sim, Val(:T_ice))
@test isapprox(Statistics.mean(T_Is), cs.model_sims.ice_sim.params.T_I_ini)
h_I = get_field(cs.model_sims.ice_sim, Val(:h_I))
@test isapprox(Statistics.mean(h_I), cs.model_sims.ice_sim.params.h_I_ini)

cs, _ = coupled_heat_equations(timestepping=:explicit)
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.u), 267.021160864323)
@test isapprox(maximum(cs.model_sims.atmos_sim.integrator.u), 270.98203813664367)
@test isapprox(minimum(cs.model_sims.atmos_sim.integrator.u), 267.0)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.u), 270.99996464720743)
@test isapprox(maximum(cs.model_sims.ocean_sim.integrator.u), 271.0)
@test isapprox(minimum(cs.model_sims.ocean_sim.integrator.u), 270.998233166346)

cs, _ = coupled_heat_equations(t_max=1000, Δt_cpl=1000, a_I=1.0, ice_model_type=:temp_feedback, timestepping=:explicit)
T_Is = get_field(cs.model_sims.ice_sim, Val(:T_ice))
@test Statistics.mean(T_Is) != cs.model_sims.ice_sim.params.T_I_ini
h_I = get_field(cs.model_sims.ice_sim, Val(:h_I))
@test isapprox(Statistics.mean(h_I), cs.model_sims.ice_sim.params.h_I_ini)

cs, _ = coupled_heat_equations(t_max=1000, Δt_cpl=1000, a_I=1.0, ice_model_type=:thickness_feedback, timestepping=:explicit)
T_Is = get_field(cs.model_sims.ice_sim, Val(:T_ice))
@test Statistics.mean(T_Is) != cs.model_sims.ice_sim.params.T_I_ini
h_I = get_field(cs.model_sims.ice_sim, Val(:h_I))
@test Statistics.mean(h_I) != cs.model_sims.ice_sim.params.h_I_ini

cs, _ = coupled_heat_equations()
@test isapprox(Statistics.mean(cs.model_sims.atmos_sim.integrator.u), 267.0215118819473)
@test isapprox(Statistics.mean(cs.model_sims.ocean_sim.integrator.u), 270.99996464720743)

cs, _ = coupled_heat_equations(a_I=1.0, ice_model_type=:thickness_feedback)
T_Is = get_field(cs.model_sims.ice_sim, Val(:T_ice))
@test Statistics.mean(T_Is) != cs.model_sims.ice_sim.params.T_I_ini
h_I = get_field(cs.model_sims.ice_sim, Val(:h_I))
@test Statistics.mean(h_I) != cs.model_sims.ice_sim.params.h_I_ini

p = SimulationParameters()
p.C_H_AO = compute_C_H_AO(p)
p.C_H_AI = compute_C_H_AI(p)
restore_physical_values!(p)
@test isapprox(compute_ϱ_ana(p), 0.003535879648086959)
