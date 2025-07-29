using QuantumOptics
using LinearAlgebra
using Plots

function TavisCummings(N, g, ω_c, Δ; spin=1//2)
    # N: number of two-level systems
    # g: coupling strength
    # ω_c: frequency of the cavity mode
    # Δ: detuning of spin modes

    # Define the Fock basis and Spin basis
    b_fock = FockBasis(1)
    b_spin = SpinBasis(spin)

    # Fundamental Operators
    a = destroy(b_fock)  # Cavity annihilation operator
    at = create(b_fock)  # Cavity creation operator
    n = number(b_fock)  # Cavity number operator

    sz = sigmaz(b_spin)  # Spin z operator
    sp = sigmap(b_spin)  # Spin raising operator
    sm = sigmam(b_spin)  # Spin lowering operator

    # Iterables
    ω_s = [ω_c - 2 * Δ + i * Δ for i = 1:N]  # Spin frequencies with detuning
    if typeof(g) != Array
        g = fill(g, N)  # If g is a single value, replicate it for all spins
    end

    # Hamiltonian
    Hspin = [ω_s[i] * sz / 2 for i = 1:N] # Spin Hamiltonian
    Hcavity = ω_c * n  # Cavity Hamiltonian
    Hint = [im*g[i]*(a⊗sp⊗I⊗I + at⊗sm⊗I⊗I + a⊗sm⊗I⊗I + at⊗sp⊗I⊗I) for i = 1:N]  # Interaction Hamiltonian
    H = Hspin[1] ⊗ Hspin[2] ⊗ Hspin[3] ⊗ Hcavity + 

    return H
end

N = 3
g = 0.1
ω_c = 1.0
Δ = 0.1

[Matrix(TavisCummings(N, g, ω_c, Δ).data)[i, i] for i = 1:4:16]