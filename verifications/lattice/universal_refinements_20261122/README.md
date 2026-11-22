# Lattice packing upper bounds in 14D, 15D and 36D

Each result is an upper bound on the occupied volume fraction for every full-rank lattice packing in the stated dimension. The table rounds upward. The replay computes the rational bounds from pinned proof inputs using exact rational arithmetic and enclosing real intervals. [Exact results](certificates.json) · [Full 9–48D comparison](../../../LATTICE_RESULTS.md)

| Dimension | Lattice density upper bound |
|---:|---:|
| 14 | 0.033996852420018639356 |
| 15 | 0.025163909222702077102 |
| 36 | 0.000027145879604620582428 |

## Reciprocal vectors in dimensions 14 and 15

Scale the lattice $L$ to shortest length one. Write $V_n$ for the volume of the unit ball and $\Delta=V_n/(2^n\operatorname{covol}(L))$. If $y$ is a shortest vector of $L^*$, it is primitive, and $L\cap y^\perp$ has covolume $|y|\operatorname{covol}(L)$. A rank-$(n-1)$ density upper bound $U_{n-1}$ therefore gives $|y|\ge 2V_{n-1}\Delta/(V_n U_{n-1})$. The bound $U_n$ in rank $n$ applies to $L^*$ with packing radius $|y|/2$. Its covolume is $1/\operatorname{covol}(L)$, which gives $|y|^2\le16(U_n\Delta/V_n^2)^{2/n}$.

For $\Delta>T$, these inequalities enclose the Fourier coordinate $s=4\pi^2|y|^2/\tau$ using $T\le\Delta\le U_n$. The rank-14 parent uses Gram determinant at most 314 in its minimum-four normalization for the Fourier case. In rank 15 the parent uses a threshold for its parity profiles. The chosen $T$ in each replay includes the corresponding parent threshold. The rank-13 and rank-14 inputs come from the pinned section proofs. [Rank-14 matching certificate](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/sections/n14/MATCHING.md) · [Native rank-15 certificate guide](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/README.md)

Write the parent's Fourier function as $e^{-s/2}q(s)$, normalized so its value at the origin is $q(0)$. In rank 14, exact Bernstein coefficients certify a contribution of at least $\mu=8\cdot10^{-7}$ from the pair $\pm y$. In rank 15, exact bounds on $e^{-s/2}q(s)$ over a rational mesh certify at least $\mu=2.6\cdot10^{-7}$ from the three pairs $\pm y,\pm2y,\pm3y$. The parent Fourier sign bounds cover all remaining vectors. Thus the Poisson summation denominator is at least $q(0)+\mu$. The resulting upper bound is $\max(T,U_n q(0)/(q(0)+\mu))$, which includes lattices of density at most $T$.

## Radial moments in dimension 36

The parent's spatial kernel uses three features $(h^6,h^9,h^{12})$ with $0\le h\le1$ on the supported nonzero lattice vectors. Let $u$ be the sum of those vectors of features. The parent completes the square in the constant spherical harmonic of its kernel, minimizing $u^{\mathsf T}Bu+2z^{\mathsf T}u$ over all real $u$. For a lattice, $u_3\le u_2$ pointwise after summation. The verifier reconstructs $B$ and $z$ and checks that the unrestricted minimizer $u_*=-B^{-1}z$ satisfies $w^{\mathsf T}u_*>0$ for $w=(0,-1,1)$. The constrained minimum is higher by

$$
\gamma=\frac{(w^{\mathsf T}u_*)^2}{w^{\mathsf T}B^{-1}w}.
$$

The parent uses the quadratic minimum to bound the spatial numerator, so it falls by $\gamma$. The parent's certified reciprocal denominator is positive and at most $q(0)$. Dividing by $q(0)$ gives a conservative deduction in the same density normalization. The replay takes the larger of the revised upper and the parent's density threshold. [Parent compact-kernel certificates](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/README.md)

## Replay

Use Python 3.11 or later with assertions enabled. From this directory, install the proof dependencies and supply the pinned `sphere-packings` checkout for 14D and the corresponding ZIP for each of the other dimensions. Keep the ZIPs intact for the archive hash checks. The run replays each parent bound and the new deduction. The underlying lower-rank exact-census inputs are documented by the parent proofs. [Pinned checkout](https://github.com/AstroNinja389/sphere-packings/tree/d7ffd0df47dfeb11814631858205f05ed73ec8ff) · [Native ZIP](https://github.com/AstroNinja389/sphere-packings/raw/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/native-proofs.zip) · [25–48D ZIP](https://github.com/AstroNinja389/sphere-packings/raw/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/n25-48-proofs-r9.zip)

```sh
python -m pip install python-flint==0.9.0 numpy scipy sympy mpmath
python -B verify.py 14 /path/to/sphere-packings
python -B verify.py 15 /path/to/native-proofs.zip
python -B verify.py 36 /path/to/n25-48-proofs-r9.zip
```
