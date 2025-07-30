using QuantumOptics
using LinearAlgebra
using Plots

function TavisCummings(N, g, ω_c, ω_s, Δ, κ, γ; spin=1//2)
    # N: number of two-level systems
    # g: coupling strength
    # ω_c: frequency of the cavity mode
    # Δ: detuning of spin modes

    # Define the Fock basis and Spin basis
    b_fock = FockBasis(1)
    global b_spin = [SpinBasis(spin), SpinBasis(spin), SpinBasis(spin)]

    # Fundamental Operators
    a = destroy(b_fock)  # Cavity annihilation operator
    at = create(b_fock)  # Cavity creation operator
    n = number(b_fock)  # Cavity number operator

    sz = [sigmaz(x) for x in b_spin]  # Spin z operator
    sp = [sigmap(x) for x in b_spin] # Spin raising operator
    sm = [sigmam(x) for x in b_spin]  # Spin lowering operator

    # Iterables
    ω_s = [ω_s - 2 * Δ + i * Δ + im * γ for i = 1:N]  # Spin frequencies with detuning
    if typeof(g) != Array
        g = fill(g, N)  # If g is a single value, replicate it for all spins
    end

    # Hamiltonian
    Hspin = []
    bases = []
    push!(bases, one(b_fock))  # Start with the Fock basis
    for i = 1:N
        push!(bases, one(b_spin[i]))
    end
    for i = 1:N
        bases_current = copy(bases)
        bases_current[i+1] = sz[i]/2
        while length(bases_current) != 1
            base_int = pop!(bases_current)
            bases_current[end] = bases_current[end] ⊗ base_int
        end
        push!(Hspin, ω_s[i] * pop!(bases_current))  # Spin Hamiltonian
    end

    bases[1] = n  # Update the Fock basis for the cavity
    bases_current = copy(bases)
    while length(bases_current) != 1
        base_int = pop!(bases_current)
        bases_current[end] = bases_current[end] ⊗ base_int
    end

    Hcavity = (ω_c + im * κ) * pop!(bases_current)  # Cavity Hamiltonian

    Hint = [im*g[1]*(a⊗sp[1]⊗one(b_spin[2])⊗one(b_spin[3])
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

    H = H + 1/2 * one(basis(H)) * (3* ω_s[2] + 0 * ω_c)  # Add a phase factor to the interaction term

    H_b_sub = SubspaceBasis(basis(H), 
        [fockstate(b_fock, 0) ⊗ spinup(b_spin[1]) ⊗ spindown(b_spin[2]) ⊗ spindown(b_spin[3]),
        fockstate(b_fock, 0) ⊗ spindown(b_spin[1]) ⊗ spinup(b_spin[2]) ⊗ spindown(b_spin[3]),
        fockstate(b_fock, 0) ⊗ spindown(b_spin[1]) ⊗ spindown(b_spin[2]) ⊗ spinup(b_spin[3]),
        fockstate(b_fock, 1) ⊗ spindown(b_spin[1]) ⊗ spindown(b_spin[2]) ⊗ spindown(b_spin[3]),
    ])  # Define the subspace basis

    P = projector(H_b_sub, basis(H))  # Project the Hamiltonian onto the subspace

    H = P * H * P'  # Projected Hamiltonian

    return H
end

N = 3
g = 1e7
ω_c = 1e9
δ = -0.4e9:1e6:0.4e9
Δ = 0.1e9
κ = g/4
γ = g/4

energies_1 = []
energies_2 = []
energies_3 = []
energies_4 = []

for δ in δ
    ω_s = ω_c + δ
    H = Matrix(TavisCummings(N, g, ω_c, ω_s, Δ, κ, γ).data)# + ω_s * I
    H = H .* [1 1 1 1; 1 1 1 1; 1 1 1 1; -1 -1 -1 1] # Ensure Hermitian
    energy = eigen(H).values
    push!(energies_1, real(energy[1])/ω_c)
    push!(energies_2, real(energy[2])/ω_c)
    push!(energies_3, real(energy[3])/ω_c)
    push!(energies_4, real(energy[4])/ω_c)
end 

plot()
plot!(δ/ω_c, energies_1, label="E1", color=:red)
plot!(δ/ω_c, energies_2, label="E2", color=:blue)
plot!(δ/ω_c, energies_3, label="E3", color=:green)
plot!(δ/ω_c, energies_4, label="E4", color=:orange) 