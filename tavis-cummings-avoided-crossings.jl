using QuantumOptics
using LinearAlgebra, ArnoldiMethod, SparseArrays
using Plots
using BenchmarkTools
using TimerOutputs

function TavisCummings(N, g, ω_c, ω_s, κ, γ; spin=1//2, open_system=false)

    if open_system == false
        γ = 0
        κ = 0
    end

    H_matrix = zeros(Complex{Float64}, N+1, N+1)
    H_matrix[1, 1] = ω_c + im * κ
    H_matrix[1, 2:end] .= im .* g
    H_matrix[2:end, 1] .= - im .* g

    if isa(ω_s, Array) == false
        ω_s = fill(ω_s, N)
    end

    if isa(γ, Array) == false
        γ = fill(γ, N)
    end

    for i = 2:N+1
        H_matrix[i, i] = ω_s[i-1] + im * γ[i-1]
    end

    b = NLevelBasis(N+1)

    H = Operator(b, b, H_matrix)

    return H
end


N = 8
g = 1e7
ω_c = 1e9
δ = -0.1e9:1e6:0.1e9
Δ = 0.2e8
κ = g/4
γ = 0

energies = [[] for _ in 1:N+1]


for δ in δ
    ω_s = [ω_c + δ + Δ * (i-N/2) for i in 1:N]
    H = Matrix(TavisCummings(N, g, ω_c, ω_s, κ, γ; open_system=false).data)
    energy = eigvals(H)
    for i = 1:N+1
        push!(energies[i], real(energy[i]) / ω_c)
    end
end

plot()
for i = 1:N+1
    plot!(δ/ω_c, energies[i], color=:red)
end
ylabel!("\$Energy (E)/ħω_c\$")
xlabel!("\$Detuning (δ)/ω_c\$")
plot!(legend=false, title="Tavis-Cummings Model Energies", show=true)