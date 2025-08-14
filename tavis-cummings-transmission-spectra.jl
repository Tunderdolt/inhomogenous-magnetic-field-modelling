#=
This script uses the Tavis-Cummings Model to evaluate the forward voltage gain 
or S₂₁ parameter of the scattering matrix for varying values of Pump/Probe and 
Spin frequencies. This is done for a two-level system only, ie the transistion 
between the |0x> and |1x> states. However, it is important to note the script
does not allow for maximally more than 1 photon in the system.

The output is a heatmap showing the values of |S₂₁|² for the 2 varying 
dimensions.

S. J. Binns
=#

# Packages =====================================================================
using QuantumOptics
using LinearAlgebra, ArnoldiMethod, SparseArrays
using Plots
using BenchmarkTools
using TimerOutputs

# Functions ====================================================================
function TavisCummings(N, g, ω_c, ω_s, κ, γ; open_system=false)

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

# Main =========================================================================
@benchmark begin
N = 100
g = 1e7
ω_c = 1e9
ω_s = 0.5e9:4e6:1.5e9
ω_p = 0.5e9:4e6:1.5e9
κ = mean(g)/4
γ = g/4

T = []
T_int = []
T_p = []

for ω_s in ω_s
    for ω_p in ω_p
        H = TavisCummings(N, g, ω_c, ω_s, κ, γ; open_system=true)
        H_int = ω_p * I - H.data
        T_calc = (
            im * κ 
            * (
                H_int[1, 1] 
                + H_int[1, 2:end]' 
                * inv(Diagonal(H_int[2:end, 2:end])) 
                * H_int[2:end, 1]
            )^-1
        )
        push!(T_int, T_calc * conj(T_calc))
    end
    push!(T, real(log10.(T_int)))
    T_int = []
end

heatmap(
    ω_s/ω_c, 
    ω_p/ω_c, 
    stack(T),
    xlabel="\$ω_s/ω_c\$", 
    ylabel="\$ω_p/ω_c\$", 
    title="Transmission Spectrum",
    color=:thermal,
    colorbar_title="Transmission",
    show=true
)
end seconds=1000  