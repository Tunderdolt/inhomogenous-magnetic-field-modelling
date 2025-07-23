using QuantumOptics
using LinearAlgebra
using Plots

function jaynesCummingsEnergies(ω_c, ω_s, g, N_cutoff, add_phase=false)
    # Define the Fock basis and Spin basis
    b_fock = FockBasis(N_cutoff)
    b_spin = SpinBasis(1//2)

    # Fundamental Operators
    a = destroy(b_fock)  # Cavity annihilation operator
    at = create(b_fock)  # Cavity creation operator
    n = number(b_fock)  # Cavity number operator

    sz = sigmaz(b_spin)  # Spin z operator
    sp = sigmap(b_spin)  # Spin raising operator
    sm = sigmam(b_spin)  # Spin lowering operator

    # Hamiltonian
    Hatom = ω_s * sz / 2  # Spin Hamiltonian
    Hcavity = ω_c * n  # Cavity Hamiltonian
    Hint = g*(a⊗sp + at⊗sm + a⊗sm + at⊗sp)  # Interaction Hamiltonian
    H = one(b_fock) ⊗ Hatom + Hcavity ⊗ one(b_spin) + Hint

    if add_phase
        # Add a phase factor to the interaction term
        H_matrix = Matrix(H.data) + 1/2 * ω_s * I
    else
        H_matrix = Matrix(H.data)
    end

    # Diagonalize the Hamiltonian
    E, V = eigen(H_matrix)
        
    return E
end

gr()

ω_c = 1.0
g = 0.1
N_cutoff = 1
ω_s = 0.8:0.01:1.2

plus_state = []
minus_state = []
energies_1 = []
energies_2 = []

for ω_s in 0.8:0.01:1.2
    energy_1 = 1/2 * (ω_s + ω_c + sqrt((ω_s - ω_c)^2 + 4 * g^2))
    energy_2 = 1/2 * (ω_s + ω_c - sqrt((ω_s - ω_c)^2 + 4 * g^2))
    push!(plus_state, energy_1)
    push!(minus_state, energy_2)
    energy = jaynesCummingsEnergies(ω_c, ω_s, g, N_cutoff, true)
    push!(energies_1, energy[2])
    push!(energies_2, energy[3])
end 

plot()
plot!(
    ω_s/ω_c, 
    energies_1, 
    label="|+>", 
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
    plus_state, 
    label="|+>", 
    color=:red,
)
plot!(
    ω_s/ω_c, 
    minus_state, 
    label="|->", 
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