# A strict determinant bound in dimension 13

For a rank-13 lattice of squared minimum four write $D$ for its Gram determinant. Its packing density is $V_{13}/\sqrt D$. The rank-12 determinant floor 729 and Mordell's inequality give $D\ge D_0=(729/4)^{13/11}$. We give a rational number $r>D_0$ with $D\ge r$ for a densest rank-13 lattice. The resulting bound on the density, rounded upward to 20 significant figures, is **0.042023513562577391279**. The reduction from the **exact** Mordell density is approximately **$8.03\times10^{-636}$ percent**. The two bounds agree at every decimal place shown. [Benchmark and its rank-12 input](https://github.com/AstroNinja389/sphere-packings/blob/7b63ee2094b7489aea3da83349c547c8da4cb388/sections/n13/SECTION-TRANSFER.md)

## A denominator bound for an optimal lattice

A densest lattice exists by compactness after fixing the minimum. It is perfect. Indeed, if its minimum-vector dyads did not span the symmetric matrices, a nonzero symmetric perturbation orthogonal to them would preserve all minimum lengths in both directions. Positive definiteness and the discrete gap to the next shell keep both sufficiently small perturbations admissible. Strict concavity of the log determinant would then give a smaller determinant in one direction.

Write $B$ for any basis matrix of this lattice, $Q=B^{\mathsf T}B$ and $D=\det Q$. Choose 91 minimum vectors $y_i=Bx_i$ whose dyads are linearly independent. Make an integer matrix $A$ from the 91 rows of monomials $x_{ij}^2$ and $x_{ij}x_{ik}$ for $j<k$. Its determinant has positive absolute value $b$. The equations $Q[x_i]=4$ give $As=4\mathbf1$, where $s$ lists the diagonal entries of $Q$ and twice its off-diagonal entries. By Cramer's rule $bQ$ is integral. Thus the reduced denominator of $D$ is at most $b^{13}$.

In orthonormal coordinates on symmetric matrices a dyad $yy^{\mathsf T}$ has diagonal coordinates $y_j^2$ and off-diagonal coordinates $\sqrt2\,y_j y_k$. The 78 factors of $\sqrt2$ and the change $X\mapsto BXB^{\mathsf T}$ multiply the absolute row determinant by $2^{39}|\det B|^{14}$. The latter exponent follows by taking singular values $\sigma_j$ of $B$ and multiplying the factors $\sigma_j\sigma_k$ for $j\le k$. Each transformed row has Frobenius norm $\|y_i\|^2=4$. Hadamard's inequality now gives

$$
b\,2^{39}D^7 \le 4^{91}.
$$

Since $D\ge D_0>469$, the integer $b$ is at most $H=\lfloor4^{91}/(2^{39}469^7)\rfloor$. Every candidate for the optimum therefore has determinant denominator at most $H^{13}$.

## An explicit strict improvement

The verifier finds neighboring fractions $a/b'<D_0<c/d$ with $cb'-ad=1$ and $b'+d>H^{13}$. A reduced fraction strictly between them has denominator at least $b'+d$. Indeed, for such a fraction $p/q$ the identity $q=d(b'p-aq)+b'(cq-dp)$ is a sum of positive integer multiples of $d$ and $b'$. Consequently the optimum has $D\ge c/d$. The checker verifies both comparisons with $D_0$ using eleventh powers in exact rational arithmetic. It then encloses $V_{13}/\sqrt{c/d}$ and the exact benchmark at sufficient precision to distinguish them.

## Replay

Run the following from this directory. The output includes the neighboring fractions, the rational density upper bound and a directed interval for the percentage improvement.

```sh
python -m pip install python-flint==0.9.0
python -B verify.py
```
