export SimulationParameters, restore_physical_values!, get_vertical_space

Base.@kwdef mutable struct SimulationParameters{FT}
    a_I::FT = 0.0
    ρ_A::FT = 1.225
    ρ_O::FT = 1e3
    c_A::FT = 1005.0
    c_O::FT = 4182.0
    ν_O::FT = 1e-6
    ν_A::FT = 1.5e-5
    k_A::FT = 0.02364
    k_O::FT = 0.58
    k_I::FT = 2.03
    alb_I::FT = 0.8
    A::FT = 309.8
    B::FT = 3.69
    ϵ::FT = 0.98
    LW_in = t -> FT(150.0)
    SW_in = t -> FT(200.0)
    J_q = t -> FT(0.0)
    J_s = nothing
    T_Ib::FT = 273 - 1.8
    q_I::FT = 2.75e8
    C_H_IO::FT = 5e-3
    C_H_AI::FT = 1.4e-3
    C_H_AO::FT = 1e-3
    Δu_AO::FT = 4.0
    Δu_AI::FT = 5.0
    Δu_IO::FT = 1.0
    h_O::FT = -50.0
    h_A::FT = 200.0
    z_A0::FT = 10.0
    z_O0::FT = 1.0
    T_A_ini::FT = 267.0
    T_O_ini::FT = 271.0
    T_I_ini::FT = 270.0
    h_I_ini::FT = 1.0
    t_0::FT = 0.0
    t_max::FT = 3600.0
    Δt_cpl::FT = 100.0
    Δt_min::FT = 1.0
    n_t_A = 50
    n_t_O = 1
    n_t_I = 1
    n_A = 200
    n_O = 50
    α_O::FT = k_O / (ρ_O * c_O)
    α_A::FT = k_A / (ρ_A * c_A)
    C_AO::FT = ρ_A * c_A * C_H_AO * Δu_AO
    C_AI::FT = ρ_A * c_A * C_H_AI * Δu_AI
    C_IO::FT = ρ_O * c_O * C_H_IO * Δu_IO
    T_A = nothing
    T_O = nothing
    F_AO = nothing
    F_O = nothing
    T_Is = nothing
    stable_range = nothing
    ice_model_type = :constant
    timestepping = :implicit
end

Base.eltype(::SimulationParameters{FT}) where {FT} = FT

function restore_physical_values!(p::SimulationParameters)
    p.α_O = p.k_O / (p.ρ_O * p.c_O)
    p.α_A = p.k_A / (p.ρ_A * p.c_A)
    p.C_AO = p.ρ_A * p.c_A * p.C_H_AO * p.Δu_AO
    p.C_AI = p.ρ_A * p.c_A * p.C_H_AI * p.Δu_AI
    p.C_IO = p.ρ_O * p.c_O * p.C_H_IO * p.Δu_IO
end

function get_vertical_space(lower_boundary::FT, upper_boundary::FT, nelems::Int) where {FT}
    device = Utilities.get_device(Dict("device" => "auto"))
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(lower_boundary),
        CC.Geometry.ZPoint{FT}(upper_boundary);
        boundary_names=(:bottom, :top),
    )
    mesh = CC.Meshes.IntervalMesh(domain, nelems=nelems)
    return CC.Spaces.CenterFiniteDifferenceSpace(device, mesh)
end