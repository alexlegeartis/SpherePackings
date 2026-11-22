# Compact-plane refits for lattice packing

The certificate applies to all full-rank lattice packings in its dimension. The quantities in the table are occupied volume fractions.

| Dimension | New upper bound | Pinned parent upper bound | Relative reduction |
|---:|---:|---:|---:|
| 38 | 0.00001229554573503757040 | 0.00001234290484423536589 | 0.38369500% |

## Offline replay

Use Python 3.11 or later, NumPy, and python-flint 0.9.0. The verifier accepts the [pinned parent proof ZIP](https://github.com/AstroNinja389/sphere-packings/raw/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/n25-48-proofs-r9.zip) or its extracted directory. It checks the source manifest and replays the parent ranks used for comparison before it checks the refit.

~~~sh
python -m pip install numpy python-flint==0.9.0
python -B verify.py /path/to/n25-48-proofs-r9.zip
~~~

The witness stores rational Laguerre coefficients, compact harmonic matrices, and a compact origin vector. The replay uses rational Bernstein subdivision for the plane signs, exact polynomial arithmetic for the Fourier signs, and Arb bounds for the exponential and gamma factors. It needs no optimizer.

## Density inequality

Normalize the squared minimum of the lattice to one. For squared length $a$ in $[1,R]$, let $H(a)$ have entries $((R-a)/(R-1))^{k_j}$, using the powers in the witness. Set $H(a)=0$ for $a>R$. The kernel on nonzero vectors is

~~~math
K(v,w)=\sum_{j\geq0}H(\lVert v\rVert^2)^{\mathsf T} A_j
  H(\lVert w\rVert^2)G_{2j}\!\left(\frac{\langle v,w\rangle}
  {\lVert v\rVert\lVert w\rVert}\right).
~~~

Each $A_j$ is positive semidefinite and each $G_{2j}$ is the normalized degree-$2j$ zonal polynomial in the dimension under consideration. The spherical addition formula makes all nonconstant harmonic contributions to the kernel sum nonnegative. For $u=\sum_{v\ne0}H(\lVert v\rVert^2)$, completing the square bounds the constant contribution with its origin terms by

~~~math
s+2z^{\mathsf T}u+u^{\mathsf T}A_0u
  \geq 0,\qquad s=z^{\mathsf T}A_0^{-1}z.
~~~

The verifier bounds the off-diagonal kernel sum plane by plane. A rank-two plane has at most three supported opposite pairs, since two distinct pairs with the same nonzero class modulo twice the plane lattice would yield two nonzero half-sums with total squared length less than two. A pair of independent supported vectors generates the plane lattice because its area is at most $R<\sqrt3$, whereas every lattice of squared minimum at least one in that plane has covolume at least $\sqrt3/2$. In a two-line plane both the sum and difference of the generating vectors lie at or outside the support boundary. In a three-line plane the third pair is represented by their sum or difference. The exact sign checks cover the corresponding two-line interaction and three-line interaction sum. Hence the off-diagonal sum is nonpositive. For each nonzero supported vector $v$, the terms at $w=v$ and $w=-v$ both contribute the diagonal value. This gives

~~~math
0\le s+\sum_{v\ne0}
   \bigl(2z^{\mathsf T}H(\lVert v\rVert^2)+2K(v,v)\bigr).
~~~

The radial certificate bounds $P(\tau a)e^{-\tau a/2}+2z^{\mathsf T}H(a)+2K(v,v)$ by zero for $a\ge1$. This bounds the sum of the primal Gaussian polynomial on the lattice by $P(0)+s$.

The reciprocal polynomial $Q$ is nonnegative for every squared length. Its Fourier sum on the reciprocal lattice is therefore at least its origin term. Poisson summation bounds the density by

~~~math
\frac{(\tau/8)^{n/2}}{\Gamma(1+n/2)}
\frac{P(0)+s}{Q(0)}.
~~~

The benchmark is the exact final bound in the pinned parent archive. [Parent bounds](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/bounds-25-48.json)
