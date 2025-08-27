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
using ProgressMeter

# Functions ====================================================================
function tavis_cummings(N, g, ω_c, ω_s, κ, γ; open_system=false)
    """Returns the Hmailtonian described by the Tavis-Cummings model for a
    system with N spins and maximally 1 photon with regards to the single
    excitation subspace.

    Parameters
    ----------
    N :: Int
        The number of spins in the system
    g :: Float
        The coupling strength between the spins and the cavity
    ω_c :: Float
        The cavity frequency (Hz)
    ω_s :: Float or Array{Float64}
        The spin frequency (Hz) or an array of frequencies for each spin
    κ :: Float
        The cavity dissipation rate (Hz)
    γ :: Float or Array{Float64}
        The spin dissipation rate (Hz) or an array of rates for each spin
    open_system :: Bool = false
        Condition for the system being open and allowing dissipation
    
    Returns
    -------
    H :: Operator
        The Hamiltonian operator for the system
    """
    if !open_system
        γ = 0
        κ = 0
    end

    # Ensure ω_s and γ are arrays
    ω_s = isa(ω_s, AbstractArray) ? ω_s : fill(ω_s, N)
    γ   = isa(γ, AbstractArray) ? γ   : fill(γ, N)
    # Collect entries as triplets
    rows = Int[]
    cols = Int[]
    vals = ComplexF64[]
    # cavity self-energy
    push!(rows, 1); push!(cols, 1); push!(vals, ω_c + im*κ)
    # cavity–spin couplings (row 1, cols 2:N+1)
    append!(rows, fill(1, N))
    append!(cols, 2:N+1)
    append!(vals, fill(im*g, N))
    # spin–cavity couplings (rows 2:N+1, col 1)
    append!(rows, 2:N+1)
    append!(cols, fill(1, N))
    append!(vals, fill(-im*g, N))
    # diagonal spin terms
    append!(rows, 2:N+1)
    append!(cols, 2:N+1)
    append!(vals, ω_s .+ im .* γ)
    # Build sparse matrix in one pass
    H_matrix = sparse(rows, cols, vals, N+1, N+1)
    # Wrap in operator
    b = NLevelBasis(N+1)
    return Operator(b, b, H_matrix)
end

# Main =========================================================================
g = 1e7
ω_c = 1e9
ω_s = -100e9:840e6:110e9
ω_p = -100e9:840e6:110e9
ω_s = 0e9:8e6:2e9
ω_p = 0e9:8e6:2e9
κ = mean(g)/4
γ = g/4

T = []
T_int = []
T_p = []


@showprogress for ω_s in ω_s
    H = tavis_cummings(Int(100), g, ω_c, ω_s, κ, γ; open_system=true)
    for ω_p in ω_p
        H_int = ω_p * I - H.data
        T_calc = (
            im * κ 
            * (
                H_int[1, 1] 
                - 
                sum(H_int[1, 2:end] 
                ./ diag(H_int[2:end, 2:end])
                .* H_int[2:end, 1])
            )^-1
        )
        push!(T_int, T_calc * conj(T_calc))
    end
    push!(T, real(log10.(T_int)))
    T_int = []
end
display(
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
)