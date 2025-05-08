function Ψ_A_T(ζ)
    if ζ > 0
        a = 1.0
        b = 2 / 3
        c = 5.0
        d = 0.35
        return -b * (ζ - (c / d)) * exp(-d * ζ) - b * c / d - (1 + b * a * ζ)^1.5 + 1
    end
    x = (1 - 16 * ζ)^(1 / 4)
    return 2 * log(0.5 * (1 + x^2))
end

function Ψ_A_u(ζ)
    if ζ > 0
        a = 1.0
        b = 2 / 3
        c = 5.0
        d = 0.35
        return -b * (ζ - (c / d)) * exp(-d * ζ) - b * c / d - a * ζ
    end
    x = (1 - 16 * ζ)^(1 / 4)
    return π / 2 - 2 * atan(x) + log((1 + x)^2 * (1 + x^2) / 8)
end

function Ψ_O_unstable(x)
    return sqrt(3) * (atan(sqrt(3)) - atan(1 / sqrt(3) * (2 * x + 1))) +
           (3 / 2) * log((x^2 + x + 1) / 3)
end

function Ψ_O_T(ζ)
    if ζ > 0
        return -5 * ζ
    end
    x = (1 - 25 * ζ)^(1 / 3)
    return ψ_O_unstable(x)
end

function Ψ_O_u(ζ)
    if ζ > 0
        return -5 * ζ
    end
    x = (1 - 14 * ζ)^(1 / 3)
    return ψ_O_unstable(x)
end

function z_ruAI(a_i::Float64)
    return max(1e-3, 0.93e-3 * (1 - a_i) + 6.05e-3 * exp(-17 * (a_i - 0.5)^2))
end

function compute_C_H_AI(p::SimulationParameters; κ=0.4, L_AI=50.0, z_rTAI=1e-3)
    momentum_contribution = log(p.z_0numA / z_ruAI(p.a_i)) - Ψ_A_u(p.z_0numA / L_AI)
    heat_contribution = log(p.z_0numA / z_rTAI) - Ψ_A_T(p.z_0numA / L_AI)
    return κ^2 / (momentum_contribution * heat_contribution)
end

function compute_C_H_AO(p::SimulationParameters; α_eos=1.8e-4, κ=0.4, L_AO=50.0, z_ruAO=2e-4, z_rTAO=2e-4)
    λ_u = sqrt(p.ρ_atm / p.ρ_oce)
    λ_T = λ_u * p.c_atm / p.c_oce
    μ = p.ν_O / p.ν_A
    L_OA = λ_u^2 / (p.T_atm_ini * α_eos * λ_T)
    momentum_atm = log(p.z_0numA / z_ruAO) - Ψ_A_u(p.z_0numA / L_AO)
    momentum_oce = log(λ_u * p.z_0numO / (z_ruAO * μ)) - Ψ_O_u(p.z_0numO / L_OA)
    momentum_contribution = momentum_atm + λ_u * momentum_oce

    heat_atm = log(p.z_0numA / z_rTAO) - Ψ_A_T(p.z_0numA / L_AO)
    heat_oce = log(λ_T * p.z_0numO / (z_rTAO * μ)) - Ψ_O_T(p.z_0numO / L_OA)
    heat_contribution = heat_atm + λ_T * heat_oce
    return κ^2 / (momentum_contribution * heat_contribution)
end
