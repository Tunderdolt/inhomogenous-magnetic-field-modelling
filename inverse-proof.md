### Prove 
```math
(\omega_p I - \mathcal{H})^{-1}_{ij=11} = \begin{pmatrix} \alpha & \mathcal{B} \\\ \Gamma & \Delta \end{pmatrix}^{-1}_{ij=11} = \left(\alpha - \sum^{N}_{i=1} \frac{\mathcal{B}_i\Gamma_i}{\Delta_{ii}}\right)^{-1}
```

### Where $\mathcal{H}$ is a Hamiltonian, and

$$\begin{pmatrix} \alpha & \mathcal{B} \\\ \Gamma & \Delta \end{pmatrix}$$ 

### is a block matrix.
##

Our Hamiltonian can be represented as:

$$\mathcal{H}=\begin{pmatrix} a & B \\\ C & D \end{pmatrix}$$

Where $a$ is a scalar, $B$ is a row vector, $C$ is a column vector, $C = B^{\dagger}$, and $D$ is a diagonal matrix. This is known as a block matrix. Please note that this is only true for our single excitation subspace.

We are strictly concerned with $inv(\omega_p\mathbf{I} - \mathcal{H})_{ij=11}$, which can still be represented as the same block matrix, with constants of the same type.

Let

$$ M = \begin{pmatrix} A_n & B \\\ C & D_m \end{pmatrix} $$

Where $A_n$ is a square matrix of dimensions $n \times n$, $B$ is a matrix of size $n \times m$, $C$ is a matrix of size $m \times n$, and $D_m$ is a square matrix of size $m \times m$.

Also let us define the upper right transformation matrix as

$$ R = \begin{pmatrix} I_n & 0 \\\ -D_{m}^{-1}C & I_m \end{pmatrix} $$

Such that

$$ MR = \begin{pmatrix} A_n & B \\\ C & D_m \end{pmatrix} \begin{pmatrix} I_n & 0 \\\ -D_{m}^{-1}C & I_m \end{pmatrix} $$

$$ = \begin{pmatrix} A_n - BD_{m}^{-1}C & B \\\ C - D_{m}D_{m}^{-1}C & D_m \end{pmatrix} $$

$$ = \begin{pmatrix} A_n - BD_{m}^{-1}C & B \\\ 0 & D_m \end{pmatrix} $$

By the distributivity of determinants

$$ det(MR) = det(M)det(R) $$

$$ det(R) = det(I_{n}I_{m}) = det(I_{n})det(I_{m}) = 1 $$

$$ \therefore \ det(M) = det(MR) $$

$$ \implies det \begin{pmatrix} A_n & B \\\ C & D_n \end{pmatrix} = det \begin{pmatrix} A_n - BD_{m}^{-1}C & B \\\ 0 & D_m \end{pmatrix} $$

$$ = det \left[ \left(A_n - BD_{m}^{-1}C \right) D_{m} \right]$$

$$ \boxed{det(M) = det \left(A_n - BD_{m}^{-1}C \right) det \left(D_{m}\right)} $$

Now we have an expression for the determinant of a block matrix, let us consider the system defined by

$$\mathcal{M} = \omega_p I - \mathcal{H} $$

Which we will define as a block matrix in the following way

$$ \mathcal{M} = \begin{pmatrix} \alpha & \mathcal{B} \\\ \Gamma & \Delta \end{pmatrix} $$

Where $\alpha$ is a scalar, $\mathcal{B}$ is a row vector of size $1 \times \nu$, $\Gamma$ is a column vector of size $\nu \times 1$ and $\Delta$ is a diagonal square matrix of size $\nu \times \nu$

We desire a solution to

$$ \mathcal{M}_{ij=11}^{-1} = \frac{1}{det(\mathcal{M})} \times det(D) $$

We can use the previous experession we found for the determinant of a block matrix as below

$$ det(\mathcal{M}) = det\left( \alpha - \mathcal{B} \Delta^{-1} \Gamma\right) det( \Delta ) $$

To simplify this expressionto

$$ \mathcal{M}_{ij=11}^{-1} = \frac{det(\Delta)}{det\left( \alpha - \mathcal{B} \Delta^{-1} \Gamma\right) det( \Delta ) } $$

$$ \mathcal{M}_{ij=11}^{-1} = \frac{1}{det\left( \alpha - \mathcal{B} \Delta^{-1} \Gamma\right)} $$

Which, due to the form of each of these matrices, can further be simplified down to 

```math
= \left(\alpha - \sum^{N}_{i=1} \frac{\mathcal{B}_i\Gamma_i}{\Delta_{ii}}\right)^{-1}
```

QED
