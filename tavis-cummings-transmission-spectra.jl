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

    # Ensure γ, ω_s, and g are arrays
    γ   = isa(γ, AbstractArray) ? γ   : fill(γ, N)
    ω_s = isa(ω_s, AbstractArray) ? ω_s : fill(ω_s, N)
    g   = isa(g, AbstractArray) ? g   : fill(g, N)

    # Collect entries as triplets
    rows = Int[]
    cols = Int[]
    vals = ComplexF64[]

    # cavity term (row 1, col 1)
    push!(rows, 1); push!(cols, 1); push!(vals, ω_c + im*κ)

    # cavity–spin couplings (row 1, cols 2:N+1)
    append!(rows, fill(1, N))
    append!(cols, 2:N+1)
    append!(vals, im.*g)

    # spin–cavity couplings (rows 2:N+1, col 1)
    append!(rows, 2:N+1)
    append!(cols, fill(1, N))
    append!(vals, -im.*g)

    # diagonal spin terms (rows 2:N+1, cols 2:N+1)
    append!(rows, 2:N+1)
    append!(cols, 2:N+1)
    append!(vals, ω_s .+ im .* γ)

    H_matrix = sparse(rows, cols, vals, N+1, N+1)

    b = NLevelBasis(N+1)

    return Operator(b, b, H_matrix)
end

function rand_q_gaussian(μ, q, std)
    """Returns one random sample from a q-Gaussian distribution with mean μ,
    entropic index q, and variance varience from a scaled Student’s 
    t-distribution.
    Parameters
    ----------
    μ :: Float
        The mean of the distribution
    q :: Float
        The entropic index of the distribution
    varience :: Float
        The variance of the distribution
    Returns
    -------
    Float
        A random sample from the q-Gaussian distribution
    """
    if std == 0
        return μ
    end
    ν = (3 - q) / (q - 1)
    if ν <= 2
        error("Degrees of freedom must be greater than 2 for finite variance.")
    end
    σ = sqrt(std^2 * (ν - 2) / ν)
    d = LocationScale(μ, σ, TDist(ν))
    return rand(d)
end

# Main =========================================================================
N = Int(2e3)
g = 1e7
ω_c = 1e9
min_freq = -1e9
max_freq = 2*ω_c - min_freq
resolution = 250
κ = mean(g)/4
γ = g/4

std = 0
μ = 0
const q = 1.39
freq_var = [rand_q_gaussian(μ, q, std) for _ in 1:N]

T = []
T_int = []
T_p = []

@showprogress for ω_s in range(min_freq, max_freq, length=resolution)
    ω_s = ω_s .+ freq_var
    H = tavis_cummings(N, g, ω_c, ω_s, κ, γ; open_system=true)
    for ω_p in range(min_freq, max_freq, length=resolution)
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
        range(min_freq, max_freq, length=resolution)/ω_c, 
        range(min_freq, max_freq, length=resolution)/ω_c, 
        stack(T),
        xlabel="\$ω_s/ω_c\$", 
        ylabel="\$ω_p/ω_c\$", 
        title="Transmission Spectrum",
        color=:thermal,
        colorbar_title="Transmission",
        show=true,
        clim=(-10, 0)
    )
)