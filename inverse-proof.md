### Claim 
```math
T \equiv (\omega_p I - \mathcal{H})^{-1}_{ij=11} \equiv \begin{pmatrix} \alpha & B \\\ \Gamma & \Delta \end{pmatrix}^{-1}_{ij=11} = \left(\alpha - \sum^{N}_{i=1} \frac{B_i\Gamma_i}{\Delta_{ii}}\right)^{-1}
```

### Where $\omega_p \ \epsilon \ &#8477;$, $\mathcal{H}$ is the Tavis-Cummings Hamiltonian in the single excitation subspace, and

$$\begin{pmatrix} \alpha & B \\\ \Gamma & \Delta \end{pmatrix}$$ 

### is a block matrix.
##
### Preliminaries

For the single excitation subspace, $\mathcal{H}$ (and hence $\omega_pI - \mathcal{H}$) has an arrowhead structure, i.e. the first row and column and the leading diagonal are the only populated elements within the system

```math
\mathcal{H} \;=\;
\begin{pmatrix}
a & b_2 & b_3 & b_4 & \cdots & b_n \\
c_2 & d_2 & 0 & 0 & \cdots & 0\\
c_3 & 0 & d_3 & 0 & \cdots & 0\\
c_4 & 0 & 0 & d_4 & \cdots & 0\\
\vdots & \vdots & \vdots & \vdots & \ddots & 0 \\
c_n & 0 & 0 & 0 & 0 & d_n
\end{pmatrix}
```

It is convenient for us to define the following block matrix

$$ \Theta = \omega_pI - \mathcal{H} = \begin{pmatrix} \alpha & B \\\ \Gamma & \Delta \end{pmatrix} $$

$\alpha$ is a scalar, $B$ is a row vector of size $1 \times N$, $\Gamma$ is a column vector or size $N \times 1$, and $\Delta$ is a square diagonal matrix of size $N \times N$, where N is the number of spins in the system.

Now, we want to find $T$, which is defined as the following

```math
T = \left(\omega_pI - \mathcal{H}\right)^{-1}_{ij=11} = \begin{pmatrix} \alpha & B \\\ \Gamma & \Delta \end{pmatrix}^{-1}_{ij=11} = \frac{det(\Delta)}{det(\Theta)}
```

Now, $det(\Delta)$ is trivial as $\Delta$ is a diagonal matrix, however finding an expression for $det(\Theta)$ is more challenging.

##
### Derivation of the Schur Complement

Let

$$ \mathcal{M} = \begin{pmatrix} \mathcal{A} & \mathcal{B} \\\ \mathcal{C} & \mathcal{D} \end{pmatrix} $$

Where $\mathcal{A}$ is a square matrix of dimensions $n \times n$, $\mathcal{B}$ is a matrix of size $n \times m$, $\mathcal{C}$ is a matrix of size $m \times n$, and $\mathcal{D}$ is a square matrix of size $m \times m$.

Also let us define the upper right transformation matrix as

$$ \mathcal{R} = \begin{pmatrix} I_n & 0 \\\ -\mathcal{D}^{-1}\mathcal{C} & I_m \end{pmatrix} $$

Such that

$$ \mathcal{MR} = \begin{pmatrix} \mathcal{A} & \mathcal{B} \\\ \mathcal{C} & \mathcal{D} \end{pmatrix} \begin{pmatrix} I_n & 0 \\\ -\mathcal{D}^{-1}\mathcal{C} & I_m \end{pmatrix} $$

$$ = \begin{pmatrix} \mathcal{A} - \mathcal{BD}^{-1}\mathcal{C} & \mathcal{B} \\\ \mathcal{C} - \mathcal{D}\mathcal{D}^{-1}\mathcal{C} & \mathcal{D} \end{pmatrix} $$

$$ = \begin{pmatrix} \mathcal{A} - \mathcal{BD}^{-1}\mathcal{C} & \mathcal{B} \\\ 0 & \mathcal{D} \end{pmatrix} $$

By the distributivity of determinants

$$ det(\mathcal{MR}) = det(\mathcal{M})det(\mathcal{R}) $$

$$ det(\mathcal{R}) = det(I_{n}I_{m}) = det(I_{n})det(I_{m}) = 1 $$

$$ \therefore \ det(\mathcal{M}) = det(\mathcal{MR}) $$

$$ \implies det \begin{pmatrix} \mathcal{A} & \mathcal{B} \\\ \mathcal{C} & \mathcal{D} \end{pmatrix} = det \begin{pmatrix} \mathcal{A} - \mathcal{BD}^{-1}\mathcal{C} & \mathcal{B} \\\ 0 & \mathcal{D} \end{pmatrix} $$

$$ = det \left[ \left(\mathcal{A} - \mathcal{BD}^{-1}\mathcal{C} \right) \mathcal{D} \right]$$

$$ \boxed{det(\mathcal{M}) = det \left(\mathcal{A} - \mathcal{BD}^{-1}\mathcal{C} \right) det \left(\mathcal{D}\right)} $$

Which is the Schur Complement.

##

Now, returning back to the previous expression

$$ T = \frac{det(\Delta)}{det(\Theta)} $$

Using the Schur Complement, it can be shown that

$$ det(\Theta) = det\left( \alpha - B \Delta^{-1} \Gamma\right) det( \Delta ) $$

To simplify this expression to

$$ T = \frac{det(\Delta)}{det\left( \alpha - B \Delta^{-1} \Gamma\right) det( \Delta ) } $$

$$ T = \frac{1}{det\left( \alpha - B \Delta^{-1} \Gamma\right)} $$

Which, due to the form of each of these matrices, can further be simplified down to 

```math
= \left(\alpha - \sum^{N}_{i=1} \frac{B_i\Gamma_i}{\Delta_{ii}}\right)^{-1}
```

QED
