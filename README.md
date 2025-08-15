# inhomogenous-magnetic-field-modelling

The goal of this project is to model the inhomogeniety of the magnetic field within a nitrogen doped diamond crystal for use in room temperature MASER technologies. 

To achieve this we have considered 2 parts: 
1. To consider the system as time-independent Hamiltonians and use this information to calculate the allowable states, or eigeneneries, as well as the variation in the S<sub>21</sub> parameter for &#969;<sub>s</sub> and &#969;<sub>p</sub>, for single and multi-spin systems using the Jaynes-Cummings and Tavis-Cummings model respecitvely.
2. To take the Hamiltonian for the system and use the Heisenberg Equation of motion to find the time evolution of the system.

Most of this has been done previously, however for most of the previous work it has been modelling a homogenous field/constant value of g, the interaction strength. It is also the case that most of this code has been written within academia and much of the source code is not publically available, so the hope is to allow others to have use of this code and be able to modify it accordingly for their needs.

## Functions
All scripts that use the same model use the same functions which will be described here.

### jaynesCummingsEnergies(&#969;<sub>c</sub>, &#969;<sub>s</sub>, g, N_cutoff, &#954;, &#611;; add_phase=false, open_system=false, return_matrix=false)
This function describes the Jaynes-Cummings Hamiltonian as below:

$$\mathcal{H} = \hbar(\omega_{c} + i\kappa)a^{\dagger}a + \frac{\hbar}{2}(\omega_{s} + i\gamma)\sigma_{\mathcal{z}} + i\hbar g(\sigma_{+} + \sigma_{-})(a + a^\dagger)$$

Where $\omega_c$ is the cavity frequency, $\omega_s$ is the spin frequency, $g$ is the interaction strength, $\hbar$ is the reduced Planck's constant, $\kappa$ and $\gamma$ are the dissipation rates in the system, $a$ is the annihilation operator, $a^\dagger$ is the creation operator, $\sigma_{\mathcal{z}}$ is the spin operator, $\sigma_{+}$ is the raising operator, and $\sigma_{-}$ is the lowering operator.

All Operators are constant and so defined within the function with use of the [QuantumOptics.jl package]. N_cutoff defines the maximum number of photons allowed in the system. The keyword arguements define the behaviour of the function, add_phase=true will add a constant phase of $\omega_s$ to the system, open_system=true will include the dissipation rathes in the system, otherwise it will set those rates to 0, and return_matrix=true will return the matrix that defines the Hamiltonian rather than the eigenenergies.

Will standardly return eigenenergies of the system.

### TavisCummings(N, g, &#969;<sub>c</sub>, &#969;<sub>s</sub>, &#954;, &#611;; open_system=false)
This function describes the Tavis-Cummings Hamiltonian as below:

$$\mathcal{H} = \hbar(\omega_{c} + i\kappa)a^{\dagger}a + \frac{\hbar}{2}\sum_{j=1}^{N}{(\omega_{s} + i\gamma)\_{j}\sigma_{j}^{\mathcal{z}}} + i\hbar \sum_{j=1}{N} g_{j}(\sigma_{j}^{+} + \sigma_{j}^{-})(a + a^\dagger)$$

Where $\omega_c$ is the cavity frequency, $\omega_s$ is the spin frequency, $g$ is the interaction strength, $\hbar$ is the reduced Planck's constant, $\kappa$ and $\gamma$ are the dissipation rates in the system, $a$ is the annihilation operator, $a^\dagger$ is the creation operator, $\sigma_{\mathcal{z}}$ is the spin operator, $\sigma_{+}$ is the raising operator, and $\sigma_{-}$ is the lowering operator.

$N$ is the number of spins in the system, and open_system=false sets the values of $\kappa$ and $\gamma$ to 0.

However, this function does not produce the full Hilbert space, but rather the single excitation subspace. This is an approximation of the system and has been made because the full size of the Hilbert space is far too large for a large number of spins, for example as a sparse matrix containing Complex{Float64} data types, the object would be of the size 10^31 bytes, which is far too large for any currently existing hardware. Reducing to the single excitation subspace makes the size of the matrix describing the Hilbert space go from $2^{N+1}\times2^{N+1}$ to $N+1\times N+1$.

This function does apply quantum operator to it's output, but that is very limited in comparison to the jaynesCummingsHamiltonian, mainly to reduce runtime for large numbers of spins.

Returns the Hamiltonian for the Tavis-Cummings Model as an operator from the [QuantumOptics.jl package].

## Scripts
All benchmarking results have been made on my own HP Envy with 16GB of RAM and an intel evo i5 processing core.

### jaynes-cummings-energy-levels
This script uses the [jaynesCummingsEnergies](#jaynesCummingsEnergies(&#969;<sub>c</sub>,-&#969;<sub>s</sub>,-g,-N_cutoff,-&#954;,-&#611;;-add_phase=false,-open_system=false,-return_matrix=false)) function to plot the eigenergies of the system for any number of photons. N_cutoff can be changed to increase the number of states plotted.
<img width="597" height="393" alt="image" src="https://github.com/user-attachments/assets/3f87d9cc-6414-4f57-9709-8e14581a5481" />
<img width="538" height="172" alt="image" src="https://github.com/user-attachments/assets/4a07e7ef-20a3-4837-847a-255f77d04d44" />

Above is shown the output for this script as well as the benchmarking results for 1000 samples, suggesting it should be expected that the code should run within 30ms for those conditions.

### jaynes-cummings-avoided-crossings
This script works almost identically to [jaynes-cummings-energy-levels](#jaynes-cummings-energy-levels), the only differences being that the add_phase arguement is now true, and we zoom in on the avoided crossing area of the graph rather than showing the whole graph.
<img width="594" height="392" alt="image" src="https://github.com/user-attachments/assets/167bdd48-bda1-4d52-b159-37502a158203" />
<img width="534" height="171" alt="image" src="https://github.com/user-attachments/assets/d507dec3-7808-416e-8b5e-30d3e79cb0c1" />

Above is shown the output for this script as well as the benchmarking results for 1000 samples, suggesting it should be expected that the code should run within 60ms for those conditions.

### jaynes-cummings-transmission-spectra
This script calculates the forward voltage gain of the system, often represented as $S_{21}$ as described by the scattering ($\mathbf{S}$) matrix as the values of spin frequency ($\omega_s$) and probe frequency ($\omega_p$) vary and plots this as a heatmap.
<img width="592" height="396" alt="image" src="https://github.com/user-attachments/assets/bc6de100-b898-423a-893c-b6f286f74c7e" />
<img width="532" height="172" alt="image" src="https://github.com/user-attachments/assets/1658446a-f920-4fbc-bba0-ba83bdc5eeac" />

### tavis-cummings-avoided-crossings
2 spins:
<img width="529" height="178" alt="image" src="https://github.com/user-attachments/assets/af9b2408-bad2-4d69-b8d3-550ff658527d" />


[QuantumOptics.jl package]: https://qojulia.org/
