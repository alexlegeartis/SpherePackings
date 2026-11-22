# Compact-plane refits for lattice packing

The certificate applies to all full-rank lattice packings in its dimension. The quantities in the table are occupied volume fractions. Benchmark decimals are rounded upward where needed. Reductions use the exact archive values.

| Dimension | New upper bound | Pinned parent upper bound | Relative reduction |
|---:|---:|---:|---:|
| 37 | 0.00001817642454038903306 | 0.000018245562973058698414114862 | 0.37893285% |
| 38 | 0.00001229554573503757040 | 0.00001234290484423536589 | 0.38369500% |
| 39 | 0.00000845291367154137428 | 0.0000084768372418329309 | 0.28222283% |
| 42 | 0.00000346603735498588186 | 0.0000034723546108616219543444225 | 0.18193003% |
| 48 | 0.00000034510710122965769 | 3.4530564242393277928417642E-7 | 0.05749723% |

## Offline replay

Use Python 3.11 or later, NumPy, and python-flint 0.9.0. The verifier accepts the [pinned parent proof ZIP](https://github.com/AstroNinja389/sphere-packings/raw/7b63ee2094b7489aea3da83349c547c8da4cb388/proofs/n25-48-proofs-r9.zip) or its extracted directory. It checks the source manifest and replays the parent ranks used for comparison before it checks the refit.

~~~sh
python -m pip install numpy python-flint==0.9.0
python -B verify.py /path/to/n25-48-proofs-r9.zip
~~~

The witness stores rational Laguerre coefficients, compact harmonic matrices, and a compact origin vector. The replay uses rational Bernstein subdivision for the plane signs, exact polynomial arithmetic for the Fourier signs, and Arb bounds for the exponential and gamma factors. It needs no optimizer.

## How the plane kernels were found

The search varied the Fourier coefficients and all the plane matrices together. In dimension 37, reducing selected positive semidefinite blocks removed the remaining angular sign violations. In dimension 39, subtracting a nonnegative even angular polynomial that vanishes at $t=1$ tightened the off-diagonal inequality while preserving its spatial diagonal.

In dimensions 37, 42 and 48, write $A_0=zz^{\mathsf T}/s+C$ with $s=z^{\mathsf T}A_0^{-1}z$. Cauchy–Schwarz gives $C\succeq0$. Keeping the first term while reducing positive rank-one components of $C$ preserves the same Schur allowance $s$. Rank-one components of the higher harmonic blocks can also be reduced. This makes it possible to repair the plane signs with a small, explicitly calculated loss on the one component that controls the origin cost. The verifier checks the resulting rational matrices and Fourier coefficients directly.

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

The benchmark is the exact final bound in the pinned parent archive. [Parent bounds](https://github.com/AstroNinja389/sphere-packings/blob/7b63ee2094b7489aea3da83349c547c8da4cb388/proofs/bounds-25-48.json)
