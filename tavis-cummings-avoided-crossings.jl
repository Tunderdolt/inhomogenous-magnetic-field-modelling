using QuantumOptics
using LinearAlgebra
using Plots
using BenchmarkTools

function TavisCummings(N, g, ω_c, ω_s, Δ, κ, γ; spin=1//2)
    # N: number of two-level systems
    # g: coupling strength
    # ω_c: frequency of the cavity mode
    # Δ: detuning of spin modes

    # Define the Fock basis and Spin basis
    b_fock = FockBasis(1)
    b_spin = [SpinBasis(spin) for i = 1:N]

    # Fundamental Operators
    a = destroy(b_fock)  # Cavity annihilation operator
    at = create(b_fock)  # Cavity creation operator
    n = number(b_fock)  # Cavity number operator

    sz = [sigmaz(x) for x in b_spin]  # Spin z operator
    sp = [sigmap(x) for x in b_spin] # Spin raising operator
    sm = [sigmam(x) for x in b_spin]  # Spin lowering operator

    # Iterables
    ω_s = [ω_s + (i-(N+1)/2) * Δ + im * γ for i = 1:N]  # Spin frequencies with detuning
    if isa(g, Array) == false
        g = fill(g, N)  # If g is a single value, replicate it for all spins
    end

    # Hamiltonian
    Hspin = []
    bases = []

    push!(bases, one(b_fock))

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

    Hint = []
    bases_annihilation = copy(bases)
    bases_annihilation[1] = a
    bases_creation = copy(bases)
    bases_creation[1] = -at
    
    for i = 1:N
        bases_current_1 = copy(bases_annihilation)
        bases_current_2 = copy(bases_creation)
        bases_current_3 = copy(bases_annihilation)
        bases_current_4 = copy(bases_creation)
        bases_current_1[i+1] = sp[i]
        bases_current_2[i+1] = sm[i]
        bases_current_3[i+1] = sm[i]
        bases_current_4[i+1] = sp[i]
        while length(bases_current_1) != 1
            base_int_1 = pop!(bases_current_1)
            bases_current_1[end] = bases_current_1[end] ⊗ base_int_1
            base_int_2 = pop!(bases_current_2)
            bases_current_2[end] = bases_current_2[end] ⊗ base_int_2
            base_int_3 = pop!(bases_current_3)
            bases_current_3[end] = bases_current_3[end] ⊗ base_int_3
            base_int_4 = pop!(bases_current_4)
            bases_current_4[end] = bases_current_4[end] ⊗ base_int_4
        end
        push!(Hint, im * g[i] * sum([pop!(bases_current_1), pop!(bases_current_2), pop!(bases_current_3), pop!(bases_current_4)]))  # Spin raising interaction Hamiltonian
    end
    
    H = sum(Hspin) + Hcavity + sum(Hint)

    H = H + 1/2 * one(basis(H)) * (N * mean(ω_s))  # Add a phase factor to the interaction term

    reduced_hilbert_space = []
    sub_bases = []
    push!(sub_bases, fockstate(b_fock, 0))
    for i = 1:N
        push!(sub_bases, spindown(b_spin[i]))
    end

    for i = 1:N
        current_sub_space = copy(sub_bases)
        current_sub_space[i+1] = spinup(b_spin[i])
        while length(current_sub_space) != 1
            int_sub_space = pop!(current_sub_space)
            current_sub_space[end] = current_sub_space[end] ⊗ int_sub_space
        end
        push!(reduced_hilbert_space, pop!(current_sub_space))  # Reduced Hilbert space
    end

    sub_bases[1] = fockstate(b_fock, 1)

    while length(sub_bases) != 1
        int_sub_bases = pop!(sub_bases)
        sub_bases[end] = sub_bases[end] ⊗ int_sub_bases
    end

    push!(reduced_hilbert_space, pop!(sub_bases))  # Add the Fock state to the reduced Hilbert space

    H_b_sub = SubspaceBasis(basis(H), reduced_hilbert_space)  # Define the subspace basis

    P = projector(H_b_sub, basis(H))  # Project the Hamiltonian onto the subspace

    H = P * H * P'  # Projected Hamiltonian

    return H
end

N = 3
g = Array([1e7, 1e7, 1e7])
ω_c = 1e9
δ = -0.4e9:1e6:0.4e9
Δ = 0.1e9
κ = g[2]/4
γ = g[2]/4

energies = [[] for _ in 1:N+1]

t1 = time()

for δ in δ
    ω_s = ω_c + δ
    H = Matrix(TavisCummings(N, g, ω_c, ω_s, Δ, κ, γ).data)
    energy = eigvals(H)
    for i = 1:N+1
        push!(energies[i], real(energy[i]) / ω_c)
    end
end 

display(time() - t1)

plot()
for i = 1:N+1
    plot!(δ/ω_c, energies[i], label="E$i", color=:red)
end

ylabel!("\$Energy (E)/ħω_c\$")
xlabel!("\$Detuning (δ)/ħω_c\$")
plot!(legend=false, title="Tavis-Cummings Model Energies", show=true)