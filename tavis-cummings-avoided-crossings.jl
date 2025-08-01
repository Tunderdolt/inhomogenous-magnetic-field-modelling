using QuantumOptics
using LinearAlgebra
using Plots
using BenchmarkTools

function TavisCummings(N; spin=1//2)
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
        push!(Hspin, pop!(bases_current))  # Spin Hamiltonian
    end

    bases[1] = n  # Update the Fock basis for the cavity
    bases_current = copy(bases)

    while length(bases_current) != 1
        base_int = pop!(bases_current)
        bases_current[end] = bases_current[end] ⊗ base_int
    end

    Hcavity = pop!(bases_current)  # Cavity Hamiltonian

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
        push!(Hint, sum([pop!(bases_current_1), pop!(bases_current_2), pop!(bases_current_3), pop!(bases_current_4)]))  # Spin raising interaction Hamiltonian
    end

    H = Hcavity + sum(Hspin) + sum(Hint)

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

    return Dict("Hcav" => Hcavity, "Hspin" => Hspin, "Hint" => Hint, "Projector" => P)
end

function HamiltonianGenerator(N, basis, g, ω_c, ω_s_central, Δ, κ, γ)
    if isa(g, Array) == false
        g = fill(g, N)
    end

    ω_s = [ω_s_central + i * Δ for i = -N/2:1:N/2]

    Hcavity = (ω_c + im * κ) * basis["Hcav"]
    Hspin = sum([(ω_s[i] + (im * γ)) .* basis["Hspin"][i] for i = 1:N])
    Hint = sum(im .* g .* basis["Hint"])

    H = basis["Projector"] * (Hcavity + Hspin + Hint + (N/2 * ω_s_central * one(Hcavity))) * basis["Projector"]'
    
    return H
end

N = 9
g = 1e7
ω_c = 1e9
δ = -0.1e9:1e6:0.1e9
Δ = 0.2e8
κ = g/4
γ = g/4

energies = [[] for _ in 1:N+1]
H_basis = TavisCummings(N)

t1 = time()

for δ in δ
    ω_s = ω_c + δ
    H = Matrix(HamiltonianGenerator(N, H_basis, g, ω_c, ω_s, Δ, κ, γ).data)
    energy = eigvals(H)
    for i = 1:N+1
        push!(energies[i], real(energy[i]) / ω_c)
    end
end 

display(time() - t1)

plot()
for i = 1:N+1
    plot!(δ/ω_c, energies[i], color=:red)
end

ylabel!("\$Energy (E)/ħω_c\$")
xlabel!("\$Detuning (δ)/ω_c\$")
plot!(legend=false, title="Tavis-Cummings Model Energies", show=true)