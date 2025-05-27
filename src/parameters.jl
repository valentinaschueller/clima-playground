export SimulationParameters, restore_physical_values!

Base.@kwdef mutable struct SimulationParameters
    a_I::Float64 = 0.0
    ρ_A::Float64 = 1.225
    ρ_O = 1e3
    ρ_I = 917.0
    c_A = 1005.0
    c_O = 4182.0
    ν_O = 1e-6
    ν_A = 1.5e-5
    k_A = 0.02364
    k_O = 0.58
    k_I = 2.03
    alb_I = 0.8
    A = 309.8
    B = 3.69
    ϵ = 0.98
    LW_in = 150.0
    SW_in = 200.0
    T_Ib = 273.15 - 1.8
    L = 3e5
    C_H_IO = 5e-3
    C_H_AI = 1.4e-3
    C_H_AO = 1e-3
    Δu_AO = 4.0
    Δu_AI = 5.0
    Δu_IO = 1.0
    h_O = 51.0
    h_A = 210.0
    z_A0 = 10.0
    z_O0 = 1.0
    T_A_ini = 267.0
    T_O_ini = 271.0
    T_I_ini = 270.0
    h_I_ini = 1.0
    t_max::Float64 = 3600.0
    Δt_cpl::Float64 = 100.0
    Δt_min = 1.0
    n_t_A = 50
    n_t_O = 1
    n_t_I = 1
    n_A = 200
    n_O = 50
    boundary_mapping = "mean"
    α_O::Float64 = k_O / (ρ_O * c_O)
    α_A::Float64 = k_A / (ρ_A * c_A)
    C_AO::Float64 = ρ_A * c_A * C_H_AO * Δu_AO
    C_AI::Float64 = ρ_A * c_A * C_H_AI * Δu_AI
    C_IO::Float64 = ρ_O * c_O * C_H_IO * Δu_IO
    T_A = nothing
    T_O = nothing
    T_Is = nothing
    stable_range = nothing
    ice_model_type = :constant
end

function restore_physical_values!(p::SimulationParameters)
    p.α_O = p.k_O / (p.ρ_O * p.c_O)
    p.α_A = p.k_A / (p.ρ_A * p.c_A)
    p.C_AO = p.ρ_A * p.c_A * p.C_H_AO * p.Δu_AO
    p.C_AI = p.ρ_A * p.c_A * p.C_H_AI * p.Δu_AI
    p.C_IO = p.ρ_O * p.c_O * p.C_H_IO * p.Δu_IO
end

