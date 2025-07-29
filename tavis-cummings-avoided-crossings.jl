using QuantumOptics
using LinearAlgebra
using Plots

function TavisCummings(N, g, ω_c, Δ; spin=1//2)
    # N: number of two-level systems
    # g: coupling strength
    # ω_c: frequency of the cavity mode
    # Δ: detuning of spin modes

    # Define the Fock basis and Spin basis
    global b_fock = FockBasis(1)
    global b_spin = [SpinBasis(spin), SpinBasis(spin), SpinBasis(spin)]

    # Fundamental Operators
    a = destroy(b_fock)  # Cavity annihilation operator
    at = create(b_fock)  # Cavity creation operator
    global n = number(b_fock)  # Cavity number operator

    sz = [sigmaz(x) for x in b_spin]  # Spin z operator
    sp = [sigmap(x) for x in b_spin] # Spin raising operator
    sm = [sigmam(x) for x in b_spin]  # Spin lowering operator

    # Iterables
    ω_c = 1
    ω_s = [ω_c - 2 * Δ + i * Δ for i = 1:N]  # Spin frequencies with detuning
    if typeof(g) != Array
        g = fill(g, N)  # If g is a single value, replicate it for all spins
    end

    # Hamiltonian
    Hspin = [ω_s[1] * one(b_fock)⊗(sz[1]/2)⊗one(b_spin[2])⊗one(b_spin[3]),
            ω_s[2] * one(b_fock)⊗one(b_spin[1])⊗(sz[2]/2)⊗one(b_spin[3]),
            ω_s[3] * one(b_fock)⊗one(b_spin[1])⊗one(b_spin[2])⊗(sz[3]/2)]  # Expand to full tensor product space

    Hcavity = ω_c * n⊗one(b_spin[1])⊗one(b_spin[2])⊗one(b_spin[3])  # Cavity Hamiltonian

    Hint = [im*g[1]*(a⊗sp[1]⊗one(b_spin[2]) ⊗one(b_spin[3])
                    + at⊗sm[1]⊗one(b_spin[2])⊗one(b_spin[3])
                    + a⊗sm[1]⊗one(b_spin[2])⊗one(b_spin[3])
                    + at⊗sp[1]⊗one(b_spin[2])⊗one(b_spin[3])),
            im*g[2]*(a⊗one(b_spin[1])⊗sp[2]⊗one(b_spin[3])
                    + at⊗one(b_spin[1])⊗sm[2]⊗one(b_spin[3])
                    + a⊗one(b_spin[1])⊗sm[2]⊗one(b_spin[3])
                    + at⊗one(b_spin[1])⊗sp[2]⊗one(b_spin[3])), 
            im*g[3]*(a⊗one(b_spin[1])⊗one(b_spin[2])⊗sp[3]
                    + at⊗one(b_spin[1])⊗one(b_spin[2])⊗sm[3]
                    + a⊗one(b_spin[1])⊗one(b_spin[2])⊗sm[3]
                    + at⊗one(b_spin[1])⊗one(b_spin[2])⊗sp[3])]  # Interaction Hamiltonian
    
    H = sum(Hspin) + Hcavity + sum(Hint)
    H = H + 1/2 * one(basis(H)) * ω_s[2]

    H_sub = SubspaceBasis(basis(H), 
        [fockstate(b_fock, 0) ⊗ spinup(b_spin[1]) ⊗ spindown(b_spin[2]) ⊗ spindown(b_spin[3]),
        fockstate(b_fock, 0) ⊗ spindown(b_spin[1]) ⊗ spinup(b_spin[2]) ⊗ spindown(b_spin[3]),
        fockstate(b_fock, 0) ⊗ spindown(b_spin[1]) ⊗ spindown(b_spin[2]) ⊗ spinup(b_spin[3]),
        fockstate(b_fock, 1) ⊗ spindown(b_spin[1]) ⊗ spindown(b_spin[2]) ⊗ spindown(b_spin[3]),
    ])  # Define the subspace basis

    P = projector(H_sub, basis(H))  # Project the Hamiltonian onto the subspace

    H = P * H * P'  # Projected Hamiltonian

    return H
end

N = 3
g = 0.1
ω_c = 1.0
Δ = 0.1


H = reverse(TavisCummings(N, g, ω_c, Δ).data)