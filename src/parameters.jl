export SimulationParameters, restore_physical_values!

Base.@kwdef mutable struct SimulationParameters
    a_I::Float64 = 0.0
    ρ_A::Float64 = 1.225
    ρ_O::Float64 = 1e3
    c_A::Float64 = 1005.0
    c_O::Float64 = 4182.0
    ν_O::Float64 = 1e-6
    ν_A::Float64 = 1.5e-5
    k_A::Float64 = 0.02364
    k_O::Float64 = 0.58
    k_I::Float64 = 2.03
    alb_I::Float64 = 0.8
    A::Float64 = 309.8
    B::Float64 = 3.69
    ϵ::Float64 = 0.98
    LW_in = t -> 150.0
    SW_in = t -> 200.0
    J_q = t -> 0.0
    J_s = nothing
    T_Ib = 273 - 1.8
    q_I::Float64 = 2.75e8
    C_H_IO::Float64 = 5e-3
    C_H_AI::Float64 = 1.4e-3
    C_H_AO::Float64 = 1e-3
    Δu_AO::Float64 = 4.0
    Δu_AI::Float64 = 5.0
    Δu_IO::Float64 = 1.0
    h_O::Float64 = -50.0
    h_A::Float64 = 200.0
    z_A0::Float64 = 10.0
    z_O0::Float64 = 1.0
    T_A_ini = 267.0
    T_O_ini = 271.0
    T_I_ini = 270.0
    h_I_ini = 1.0
    t_0::Float64 = 0.0
    t_max::Float64 = 3600.0
    Δt_cpl::Float64 = 100.0
    Δt_min::Float64 = 1.0
    n_t_A = 50
    n_t_O = 1
    n_t_I = 1
    n_A = 200
    n_O = 50
    α_O::Float64 = k_O / (ρ_O * c_O)
    α_A::Float64 = k_A / (ρ_A * c_A)
    C_AO::Float64 = ρ_A * c_A * C_H_AO * Δu_AO
    C_AI::Float64 = ρ_A * c_A * C_H_AI * Δu_AI
    C_IO::Float64 = ρ_O * c_O * C_H_IO * Δu_IO
    T_A = nothing
    T_O = nothing
    F_AO = nothing
    F_O = nothing
    T_Is = nothing
    stable_range = nothing
    ice_model_type = :constant
    timestepping = :implicit
end

function restore_physical_values!(p::SimulationParameters)
    p.α_O = p.k_O / (p.ρ_O * p.c_O)
    p.α_A = p.k_A / (p.ρ_A * p.c_A)
    p.C_AO = p.ρ_A * p.c_A * p.C_H_AO * p.Δu_AO
    p.C_AI = p.ρ_A * p.c_A * p.C_H_AI * p.Δu_AI
    p.C_IO = p.ρ_O * p.c_O * p.C_H_IO * p.Δu_IO
end

