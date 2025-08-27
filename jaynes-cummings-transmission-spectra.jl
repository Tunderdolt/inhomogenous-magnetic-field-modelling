#=
This script uses the Jaynes-Cummings Model to evaluate the forward voltage gain 
or S₂₁ parameter of the scattering matrix for varying values of Pump and Probe
frequencies. This is done for a two-level system only, ie the transistion 
between the |0↑> and |1↓> states. However, it is important to note the script
allows for maximally more than 1 photon in the system - this can be achieved by
changing the value of N_cutoff. 

The output is a heatmap showing the values of |S₂₁|² for the 2 varying 
dimensions.

S. J. Binns
=#

# Packages =====================================================================
using QuantumOptics
using LinearAlgebra
using Plots
using BenchmarkTools

# Functions ====================================================================
function jaynes_cummings(ω_c, ω_s, g, N_cutoff, κ, γ; open_system=false)
    """Returns the Hamiltonian that describes a 2 level system for the jump of
    states between |0↑> and |1↓>.

    Parameters
    ----------
    ω_c :: Float
        Cavity/photon frequency (Hz)
    ω_s :: Float
        Spin/atom frequency (Hz)
    g :: Float
        Interaction factor
    N_cutoff :: Int
        The maximum number of photons in the system
    κ :: Float
        Cavity dissipation rate
    γ :: Float
        Spin dissipation rate
    open_system :: Bool = false
        Condition for the system being open and allowing dissipation
    
    Returns
    -------
    H :: Matrix{Complex{Float64}}
        The Hamiltonian matrix for the system
    """
    b_fock = FockBasis(N_cutoff)
    b_spin = SpinBasis(1//2)

    # Annihilation operator (a)
    a = destroy(b_fock)
    # Creation operator (a†)
    at = create(b_fock) 
    # Number operator (a†a)
    n = number(b_fock) 

    # Spin z operator (σz)
    sz = sigmaz(b_spin) 
    # Spin raising operator (σ₊)
    sp = sigmap(b_spin) 
    # Spin lowering operator (σ₋)
    sm = sigmam(b_spin)

    if open_system == false
        γ = 0
        κ = 0
    end
    
    Hatom = (ω_s + im*γ) * sz / 2 
    Hcavity = (ω_c + im*κ) * n 
    Hint = im*g*(a⊗sp - at⊗sm + a⊗sm - at⊗sp)

    H = Matrix((one(b_fock) ⊗ Hatom + Hcavity ⊗ one(b_spin) + Hint).data)

    return H
end

# Main =========================================================================
@benchmark begin
ω_c = 1.0e9
ω_s = 0.98e9:0.0001e9:1.02e9
ω_p = 0.98e9:0.0001e9:1.02e9
g = 1e6
N_cutoff = 1
κ = g/4
γ = g/4

T = []
T_int = []
T_p = []

for ω_s in ω_s
    for ω_p in ω_p
        H = jaynes_cummings(ω_c, ω_s, g, N_cutoff, κ, γ; open_system=true)
        H += 1/2 * ω_s * I
        T_calc = im * κ * inv(ω_p * I - H[4:-3:1, 4:-3:1])[1, 1]
        push!(T_int, T_calc * conj(T_calc))
    end
    push!(T, real(log10.(T_int)))
    T_int =[]
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
)
end samples=1000 seconds=2250