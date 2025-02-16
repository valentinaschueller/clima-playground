using Plots

function psi(xi, atm, stable, heat)
    if atm && stable && heat
        return -2 / 3 * (xi - (5 / 0.35)) * exp(-0.35 * xi) - (1 + (2 * xi / 3))^1.5 - (10 / 1.05) + 1
    elseif atm && stable && !heat
        return -2 / 3 * (xi - (5 / 0.35)) * exp(-0.35 * xi) - xi - (10 / 1.05)
    elseif atm && !stable && heat
        return 2 * log((1 + (1 - 16 * xi)^(1 / 2)) / 2)
    elseif atm && !stable && !heat
        x = (1 - 16 * xi)^(1 / 4)
        return pi / 2 - 2 * atan(x) + log((1 + x)^2 * (1 + x^2) / 8)
    elseif !atm && stable
        return 1 + 5 * xi
    elseif !atm && !stable
        if heat
            x = (1 - 25 * xi)^(1 / 3)
        else
            x = (1 - 14 * xi)^(1 / 3)
        end
        return sqrt(3) * (atan(sqrt(3)) - atan(1 / sqrt(3)) * (2 * x + 1)) + (3 / 2) * log((x^2 + x + 1) / 3)
    end
end

ρ_atm = Float64(1)
ρ_oce = Float64(1000)
c_atm = Float64(1000)
c_oce = Float64(4180)
lambda_u = sqrt(ρ_atm / ρ_oce)
lambda_T = sqrt(ρ_atm / ρ_oce) * c_atm / c_oce
nu_O = Float64(1e-6)
nu_A = Float64(1.5 * 1e-5)
mu = nu_O / nu_A
kappa = Float64(0.4)

z_0numA = Float64(10) # numerical zero height atmosphere [m]
z_0numO = Float64(1) # numerical zero height atmosphere [m]
z_ruAO = Float64(2 * 10^-4)
z_rTAO = Float64(2 * 10^-4)
z_rTAI = Float64(10^-3)

a_i = 0
# a_is = 0:0.01:1

# L_AO = 50
L_AOs = 10:1:200
# L_AOs = vcat(-200:1:-50, 10:1:200)

L_AI = 50
# L_AIs = 10:1:200
# L_AIs = vcat(-200:1:-50, 10:1:200) # Checking what the difference would be if we allowed for unstable conditions

# L_OA = 100
L_OAs = 30:1:550
# L_OAs = vcat(-550:1:-150, 30:1:550)
C_AOs = zeros(length(L_AOs), length(L_OAs))
C_AIs = zeros(length(a_is), length(L_AIs))
# for (i, a_i) in enumerate(a_is)
# for (j, L_AI) in enumerate(L_AIs)
for (i, L_AO) in enumerate(L_AOs)
    for (j, L_OA) in enumerate(L_OAs)
        z_ruAI = Float64(maximum([10^-3, 0.93 * 10^-3 * (1 - a_i) + 6.05 * 10^-3 * exp(-17 * (a_i - 0.5)^2)]))
        # println(z_ruAI)
        if L_AO > 0
            stable_atm = true
        else
            stable_atm = false
        end
        if L_OA > 0
            stable_oce = true
        else
            stable_oce = false
        end
        C_AO = kappa^2 / ((log(z_0numA / z_ruAO) - psi(z_0numA / L_AO, true, stable_atm, false) + lambda_u * (log(lambda_u * z_0numO / (z_ruAO * mu)) - psi(z_0numO / L_OA, false, stable_oce, false))) * (log(z_0numA / z_rTAO) - psi(z_0numA / L_AO, true, stable_atm, true) + lambda_T * (log(lambda_T * z_0numO / (z_rTAO * mu)) - psi(z_0numO / L_OA, false, stable_oce, true))))
        # C_AI = kappa^2 / (log(z_0numA / z_ruAI) - psi(z_0numA / L_AI, true, stable_atm, false)) * (log(z_0numA / z_rTAI) - psi(z_0numA / L_AI, true, stable_atm, true))
        C_AOs[i, j] = C_AO
        # C_AIs[i, j] = C_AI
    end
end

println(size(C_AOs))
# println(size(C_AIs))
# plotly()
# surface(L_AOs, L_OAs, C_AOs', xlabel="Lᴬᴼ", ylabel="Lᴼᴬ", zlabel="Cᴬᴼ")
# surface(a_is, L_AIs, C_AIs', xlabel="aᴵ", ylabel="Lᴬᴵ", zlabel="Cᴬᴵ")

println(maximum(C_AOs))
println(minimum(C_AOs))
# println(maximum(C_AIs))
# println(minimum(C_AIs))

# xlims!(10, 200)
# ylims!(30, 550)
# xlims!(0, 1)
# ylims!(10, 200)
# display(current())
# println(C_AOs)
# println(C_AO)