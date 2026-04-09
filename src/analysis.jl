export compute_ϱ_ana

function χ_A(p::SimulationParameters, s)
    return tanh(p.h_A * sqrt(s / p.α_A))
end

function χ_O(p::SimulationParameters, s)
    return tanh(p.h_O * sqrt(s / p.α_O))
end

function compute_ϱ_ana(p::SimulationParameters; s=nothing)
    if isnothing(s)
        s = im * π / p.t_max
    end
    return abs(p.k_A / p.k_O * sqrt(p.α_O / p.α_A) * χ_A(p, s) / χ_O(p, s))
end

