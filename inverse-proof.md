# Simplification of $\mathcal{H}^{-1}_{ij=11}$

Our Hamiltonian can be represented as:

$$\mathcal{H}=\begin{bmatrix} A & B \\\ C & D \end{bmatrix}$$

Where $A$ is a scalar, $B$ is a row vector, $C$ is a column vector, $C = B^{\dagger}$, and $D$ is a diagonal matrix. This is known as a block matrix. Please note that this is only true for our single excitation subspace.

We are strictly concerned with $inv(\omega_p\mathbf{I} - \mathcal{H})_{ij=11}$, which can still be represented as the same block matrix, with constants of the same type.

$$\begin{bmatrix} A & B \\\ C & D \end{bmatrix}^{-1}_{ij=11} = \frac{1}{det(AD - BC)}det(D)$$

As the determinants of matrices are commutative, and $det(M^{-1}) = det(M)^{-1}$,

$$\begin{vmatrix} A & B \\\ C & D \end{vmatrix} = det(D) \begin{vmatrix} A & B \\\ C & D \end{vmatrix} det(D^{-1})$$

$$= det(D) \begin{vmatrix} A & B \\\ C & D \end{vmatrix} \begin{vmatrix} 1 & 0 \\\ -D^{-1}C & D^{-1} \end{vmatrix}$$

$$= det(D) \begin{vmatrix} A - BD^{-1}C & BD^{-1} \\\ 0 & DD^{-1} \end{vmatrix}$$

$$= det(D)det(A - BD^{-1}C)$$

$$\therefore \\ \frac{1}{det(AD - BC)}det(D) = \frac{det(D)}{det(D)det(A - BD^{-1}C)}$$

$$= \frac{1}{det(A - BD^{-1}C)}$$

$$= \frac{1}{det(A - BD^{-1}B^\dagger)}$$

$$= (A - BD^{-1}B^\dagger)^{-1}$$

Furthermore, because D is a diagonal matrix, this can be further simplified to

```math
= A - \sum^{N}_{i=1} \frac{B_iB_i^{\dagger}}{D_{ii}}
```

This is a highly programming efiicient version of the equation especially when compared to solutions involving an inverse operation.
