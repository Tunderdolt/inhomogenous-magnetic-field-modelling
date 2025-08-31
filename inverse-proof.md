# Simplification of $\mathcal{H}^{-1}_{ij=11}$

Our Hamiltonian can be represented as:

$$\mathcal{H}=\begin{pmatrix} A & B \\\ C & D \end{pmatrix}$$

Where $A$ is a scalar, $B$ is a row vector, $C$ is a column vector, $C = B^{\dagger}$, and $D$ is a diagonal matrix. This is known as a block matrix. Please note that this is only true for our single excitation subspace.

We are strictly concerned with $inv(\omega_p\mathbf{I} - \mathcal{H})_{ij=11}$, which can still be represented as the same block matrix, with constants of the same type.

$$\begin{pmatrix} A & B \\\ C & D \end{pmatrix}^{-1}_{ij=11} = \frac{1}{det(AD - BC)}det(D)$$

As the determinants of matrices are commutative, and $det(M^{-1}) = det(M)^{-1}$,

$$det\begin{pmatrix} A & B \\\ C & D \end{pmatrix} = det(D) det\begin{pmatrix} A & B \\\ C & D \end{pmatrix} det(D^{-1})$$

$$= det(D) det\begin{pmatrix} A & B \\\ C & D \end{pmatrix} det\begin{pmatrix} 1 & 0 \\\ -D^{-1}C & D^{-1} \end{pmatrix}$$

$$= det(D) det\begin{pmatrix} A - BD^{-1}C & BD^{-1} \\\ 0 & DD^{-1} \end{pmatrix}$$

$$= det(D)det(A - BD^{-1}C)$$

$$\therefore \\ \frac{1}{det(AD - BC)}det(D) = \frac{det(D)}{det(D)det(A - BD^{-1}C)}$$

$$= \frac{1}{det(A - BD^{-1}C)}$$

$$= \frac{1}{det(A - BD^{-1}B^\dagger)}$$

$$= (A - BD^{-1}B^\dagger)^{-1}$$

Furthermore, because D is a diagonal matrix, this can be further simplified to

```math
= \left(A - \sum^{N}_{i=1} \frac{B_iB_i^{\dagger}}{D_{ii}}\right)^{-1}
```

This is a highly programming efiicient version of the equation especially when compared to solutions involving an inverse operation.

A suspicious result may be the following:

$$ det(D^{-1}) = det\begin{pmatrix} 1 & 0 \\\ -D^{-1}C & D^{-1} \end{pmatrix}$$

In order to do this we want a transformation matrix that takes our block matrix and converts it to an upper right matrix

$$ R = \begin{pmatrix} I & 0 \\\ D^{-1}C & I \end{pmatrix} $$

$$ M = det\begin{pmatrix} 1 & 0 \\\ -D^{-1}C & D^{-1} \end{pmatrix} $$

$$ MR = \begin{pmatrix} 1 & 0 \\\ -D^{-1}C & D^{-1} \end{pmatrix} \ \begin{pmatrix} I & 0 \\\ D^{-1}C & I \end{pmatrix} $$

$$ = \begin{pmatrix} I & 0 \\\ -D^{-1}C + D^{-1}C & D^{-1} \end{pmatrix} $$ 

$$ = \begin{pmatrix} I & 0 \\\ 0 & D^{-1} \end{pmatrix} $$

$$ det(MR) = det(D^{-1}) $$

But, also by the rules of multiplacativity for determinants

$$ det(MR) = det(M)det(R) $$

$$ det(R) = det(II) = 1 $$

$$ => det(MR) = det(M) $$

$$ \therefore \ det\begin{pmatrix} 1 & 0 \\\ -D^{-1}C & D^{-1} \end{pmatrix} = det(D^{-1}) $$

QED
