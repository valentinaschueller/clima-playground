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


"""Creates a dict with values and names for all variables."""
function get_var_dict()
    variable_dict = Dict(
        :Δt_cpl =>
            [Float64.([10, 100, 1000, 10000, 100000, 1000000]), L"$\Delta t_{cpl}$"],
        :h_atm => [Float64.([15, 50, 100, 200, 300, 500, 700, 1500]), L"$H^A$"],
        :h_oce => [Float64.([6, 10, 50, 100, 125, 150, 175, 200]), L"$H^O$"],
        :n_atm => [[20, 200, 2000, 20000, 200000], L"$\Delta z^A$"],
        :n_oce => [[5, 50, 500, 5000, 50000], L"$\Delta z^O$"],
        :C_H_AI => [
            Float64.(LinRange(0.0008691059360985882, 0.0027815836784748733, 10)),
            L"$C^A_I$",
        ],
        :C_H_AO => [
            Float64.(LinRange(0.0006765386884900067, 0.0014898726741724184, 10)),
            L"$C^A_O$",
        ],
        :C_H_IO => [Float64.(LinRange(0.001, 0.01, 10)), L"$C^O_I$"],
        :Δu_AO => [Float64.([0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]), L"$\Delta u$"],
        :T_atm_ini =>
            [Float64.([260, 262, 264, 266, 268, 270, 273, 276]), L"$T^A(0,z)$"],
        :T_oce_ini => [Float64.([270, 271, 272, 273, 274, 275, 276]), L"$T^O(0,z)$"],
        :T_ice_ini => [Float64.([260, 262, 264, 266, 268, 270, 273]), L"$T^I$"],
        :n_t_atm =>
            [Float64.(10 .^ LinRange(log10(1), log10(1000), 10)), L"$\Delta t^A$"],
        :n_t_oce =>
            [Float64.(10 .^ LinRange(log10(1), log10(1000), 10)), L"$\Delta t^O$"],
    )
    return variable_dict
end

"""Creates a dict with colors for the different sea ice concentrations."""
function get_color_dict()
    color_dict = Dict(
        0.0 => RGB(0.850, 0.325, 0.098),
        0.1 => RGB(0.0, 0.447, 0.741),
        0.2 => RGB(0.929, 0.694, 0.125),
        0.3 => RGB(0.494, 0.184, 0.556),
        0.4 => RGB(0.75, 0.0, 0.75),
        0.5 => RGB(0.301, 0.745, 0.933),
        0.6 => RGB(0.5, 0.5, 0.5),
        0.7 => RGB(0.466, 0.674, 0.188),
        0.8 => RGB(1.0, 0.843, 0.0),
        0.9 => RGB(0.0, 0.75, 0.75),
        1.0 => RGB(0.635, 0.078, 0.184),
    )
    linestyle_dict = Dict(
        0.0 => :solid,
        0.1 => :solid,
        0.2 => :dashdot,
        0.3 => :dash,
        0.4 => :dash,
        0.5 => :dashdotdot,
        0.6 => :solid,
        0.7 => :dot,
        0.8 => :dashdot,
        0.9 => :dashdotdot,
        1.0 => :dot,
    )
    return color_dict, linestyle_dict
end
