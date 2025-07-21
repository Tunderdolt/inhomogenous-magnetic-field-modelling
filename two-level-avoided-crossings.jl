using QuantumOptics
using LinearAlgebra
using Plots

gr()

ω_c = 1.0
g = 0.1
N_cutoff = 1
ω_s = 0.8:0.01:1.2

plus_state = []
minus_state = []

for ω_s in 0.8:0.01:1.2
    energy_1 = 1/2 * (ω_s + ω_c + sqrt((ω_s - ω_c)^2 + 4 * g^2))
    energy_2 = 1/2 * (ω_s + ω_c - sqrt((ω_s - ω_c)^2 + 4 * g^2))
    push!(plus_state, energy_1)
    push!(minus_state, energy_2)
end 

plot()
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
ylabel!("Energy (E)")
xlabel!("ω_s/ω_c", show=true, legend=false)