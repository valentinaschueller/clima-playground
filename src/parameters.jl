Base.@kwdef mutable struct SimulationParameters
    a_i::Float64 = 0.0
    ρ_atm::Float64 = 1.225
    ρ_oce = 1e3
    c_atm = 1005.0
    c_oce = 4182.0
    ν_O = 1e-6
    ν_A = 1.5e-5
    k_atm = 0.02364
    k_oce = 0.58
    C_H_IO = 5e-3
    C_H_AI = 1.4e-3
    C_H_AO = 1e-3
    Δu_AO = 4.0
    Δu_AI = 5.0
    Δu_IO = 1.0
    h_oce = 51.0
    h_atm = 210.0
    z_0numA = 10.0
    z_0numO = 1.0
    T_atm_ini = 267.0
    T_oce_ini = 271.0
    T_ice_ini = 270.0
    t_max = 3600.0
    Δt_cpl = 100.0
    Δt_min = 1.0
    n_t_atm = 50
    n_t_oce = 1
    n_atm = 200
    n_oce = 50
    boundary_mapping = "mean"
    sin_field_atm = false
    sin_field_oce = false
    α_o::Float64 = k_oce / (ρ_oce * c_oce)
    α_a::Float64 = k_atm / (ρ_atm * c_atm)
    C_AO::Float64 = ρ_atm * c_atm * C_H_AO * Δu_AO
    C_AI::Float64 = ρ_atm * c_atm * C_H_AI * Δu_AI
    C_IO::Float64 = ρ_oce * c_oce * C_H_IO * Δu_IO
end

function restore_physical_values!(p::SimulationParameters)
    p.α_o = p.k_oce / (p.ρ_oce * p.c_oce)
    p.α_a = p.k_atm / (p.ρ_atm * p.c_atm)
    p.C_AO = p.ρ_atm * p.c_atm * p.C_H_AO * p.Δu_AO
    p.C_AI = p.ρ_atm * p.c_atm * p.C_H_AI * p.Δu_AI
    p.C_IO = p.ρ_oce * p.c_oce * p.C_H_IO * p.Δu_IO
end

