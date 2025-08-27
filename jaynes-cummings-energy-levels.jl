#=
This script uses the Jaynes-Cummings Model to evaluate the eigenenergies of the 
system for varying detunings (δ = (ω_s - ω_c) / ω_c)for maximally 1 spin in the 
system and any number of photons.

The output is a graph showing the eigenenergies as a function of the detuning, 
δ.

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

# Main =========================================================================c 
ω_c = 1.0
g = 0.1
N_cutoff = 2
δ = -3:0.01:1
κ = 0
γ = 0

energies = [[] for _ in 1:(N_cutoff+1)*2]

# Flipped allows distinction between |0↑> and |1↓> states after loss of order
# from eigen decompositon
flipped = true

for δ = δ
    ω_s = 1.0 + δ
    H = jaynes_cummings(ω_c, ω_s, g, N_cutoff, κ, γ)
    energy = eigvals(H)
    if energy[1] == energy[2]
        flipped = false
    end
    if flipped
        for i = 1:2:(N_cutoff+1)*2
            push!(energies[i], energy[i+1])
            push!(energies[i+1], energy[i])
        end
    else
        for i = 1:(N_cutoff+1)*2
            push!(energies[i], energy[i])
        end
    end
end 

plot()
for i = 1:(N_cutoff+1)*2
    plot!(
        δ, 
        real(energies[i]), label= i % 2 == 1 ? "|$(Int(i/2 + 0.5))↑>" : "|$(Int(i/2))↓>", 
        color=:red, 
        linestyle= i % 2 == 1 ? :solid : :dash,
    )
end
ylabel!("\$Energy (E)/ħω_c\$")
xlabel!("\$Detuning (δ = (ω_s - ω_c) / ω_c )\$", show=true, legend=false)