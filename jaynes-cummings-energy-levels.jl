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
N_cutoff = 4
δ = -3:0.001:1

energies = [[] for _ in 1:(N_cutoff+1)*2]
println(energies)

flipped = true

for δ = δ
    ω_s = 1.0 + δ
    energy = jaynesCummingsEnergies(ω_c, ω_s, g, N_cutoff)
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
        energies[i], label= i % 2 == 1 ? "|$(Int(i/2 + 0.5))↑>" : "|$(Int(i/2))↓>", 
        color=:red, 
        linestyle= i % 2 == 1 ? :solid : :dash,
    )
end
ylabel!("\$Energy (E)/ħω_c\$")
xlabel!("\$Detuning (δ = (ω_s - ω_c) / ω_c )\$", show=true, legend=false)
