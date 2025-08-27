#=
This script uses the Jaynes-Cummings Model to evaluate the eigenenergies of the 
system for varying detunings (δ = (ω_s - ω_c) / ω_c)for maximally 1 spin in the 
system and any number of photons.

The output is a graph zoomed in on the eigenenergies around the centre of the
avoided crossing as a funciton of ω_s/ω_c
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
ω_c = 1.0
g = 0.1
N_cutoff = 1
ω_s = 0.8:0.01:1.2
κ = 0
γ = 0

energies_1 = []
energies_2 = []

for ω_s in 0.8:0.01:1.2
    H = jaynes_cummings(ω_c, ω_s, g, N_cutoff, κ, γ) + 1/2 * ω_s * I
    energy = eigvals(H)
    push!(energies_1, energy[2])
    push!(energies_2, energy[3])
end 

plot()
plot!(
    ω_s/ω_c, 
    energies_1, 
    label="|->", 
    color=:red,
)
plot!(
    ω_s/ω_c, 
    energies_2, 
    label="|+>", 
    color=:red,
)
plot!(
    ω_s/ω_c, 
    ω_s, 
    color=:black, 
    linestyle=:dash,
)
plot!(
    ω_s/ω_c, 
    [ω_c for _ in ω_s], 
    color=:black,
    linestyle=:dash,
)
ylabel!("\$Energy (E)/ħω_c\$", show=true, legend=false)
xlabel!("\$ω_s/ω_c\$", show=true, legend=false)