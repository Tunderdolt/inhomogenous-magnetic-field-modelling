using QuantumOptics
using LinearAlgebra
using Plots

function jaynesCummingsEnergies(ω_c, ω_s, g, N_cutoff, κ, γ; add_phase=false, open_system=false, return_matrix=false)
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
    Hint = g*(a⊗sp + at⊗sm + a⊗sm + at⊗sp)  # Interaction Hamiltonian
    if open_system
        # Open system Hamiltonian
        Hatom = (ω_s + im*γ) * sz / 2  # Spin Hamiltonian
        Hcavity = (ω_c + im*κ) * n  # Cavity Hamiltonian
        Hint = g*(a⊗sp + at⊗sm + a⊗sm + at⊗sp)  # Interaction Hamiltonian
        H = one(b_fock) ⊗ Hatom + Hcavity ⊗ one(b_spin) + Hint
    else
        Hatom = ω_s * sz / 2  # Spin Hamiltonian
        Hcavity = ω_c * n  # Cavity Hamiltonian
        Hint = g*(a⊗sp + at⊗sm + a⊗sm + at⊗sp)  # Interaction Hamiltonian
        H = one(b_fock) ⊗ Hatom + Hcavity ⊗ one(b_spin) + Hint
    end

    if add_phase
        # Add a phase factor to the interaction term
        H_matrix = Matrix(H.data) + 1/2 * ω_s * I
    else
        H_matrix = Matrix(H.data)
    end

    # Diagonalize the Hamiltonian
    E, V = eigen(H_matrix)
    
    if return_matrix
        return H_matrix
    else
        return E
    end
end

ω_c = 1.0
ω_s = 0.98:0.00005:1.02
ω_p = 0.98:0.00005:1.02
g = 0.002
N_cutoff = 1
κ = g/4
γ = g/4

T = []
T_int = []
T_p = []

for ω_s in ω_s
    for ω_p in ω_p
        #H = jaynesCummingsEnergies(ω_c, ω_s, g, N_cutoff, κ, γ; add_phase=true, open_system=true, return_matrix=true)
        #T_calc = tr(inv(ω_p * I - H))
        T_calc = im*κ * ((ω_p - ω_s - im*γ) / (-g^2 + (ω_p - ω_s - im*γ) * (ω_c - ω_p + im*κ)))
        push!(T_int, T_calc * conj(T_calc))
    end
    push!(T, real(log10.(T_int)))
    T_int =[]
end

heatmap(
    ω_s, 
    ω_p, 
    stack(T),
    xlabel="ω_s", 
    ylabel="ω_p", 
    title="Transmission Spectrum",
    color=:thermal,
    colorbar_title="Transmission",
)

plot!(ω_s, ω_s)