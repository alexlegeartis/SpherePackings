# Lattice packing upper bounds in 14D, 15D, 16D, 19D, 23D and 36D

Each result is an upper bound on the occupied volume fraction for every full-rank lattice packing in the stated dimension. The table rounds upward. The replay computes the rational bounds from pinned proof inputs using exact rational arithmetic and enclosing real intervals. [Exact results](certificates.json) · [Full 9–48D comparison](../../../LATTICE_RESULTS.md)

| Dimension | Lattice density upper bound |
|---:|---:|
| 14 | 0.033996852420018639356 |
| 15 | 0.025163909222702077102 |
| 16 | 0.017231432341539144664 |
| 19 | 0.0076312827957442716237 |
| 23 | 0.0025265652019047523943 |
| 36 | 0.000027145852648811524495 |

## Reciprocal vectors in dimensions 14–16

Scale the lattice $L$ to shortest length one. Write $V_n$ for the volume of the unit ball and $\Delta=V_n/(2^n\operatorname{covol}(L))$. If $y$ is a shortest vector of $L^*$, it is primitive, and $L\cap y^\perp$ has covolume $|y|\operatorname{covol}(L)$. A rank-$(n-1)$ density upper bound $U_{n-1}$ therefore gives $|y|\ge 2V_{n-1}\Delta/(V_n U_{n-1})$. The bound $U_n$ in rank $n$ applies to $L^*$ with packing radius $|y|/2$. Its covolume is $1/\operatorname{covol}(L)$, which gives $|y|^2\le16(U_n\Delta/V_n^2)^{2/n}$.

For $\Delta>T$, these inequalities enclose the Fourier coordinate $s=4\pi^2|y|^2/\tau$ using $T\le\Delta\le U_n$. The rank-14 parent uses Gram determinant at most 314 in its minimum-four normalization for the Fourier case. In rank 15 the parent uses a threshold for its parity profiles. The chosen $T$ in each replay includes the corresponding parent threshold. The rank-13 and rank-14 inputs come from the pinned section proofs. [Rank-14 matching certificate](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/sections/n14/MATCHING.md) · [Native rank-15 certificate guide](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/README.md)

Write the parent's Fourier function as $e^{-s/2}q(s)$, normalized so its value at the origin is $q(0)$. In rank 14, exact Bernstein coefficients certify a contribution of at least $\mu=8\cdot10^{-7}$ from the pair $\pm y$. In rank 15, exact bounds on $e^{-s/2}q(s)$ over a rational mesh certify at least $\mu=2.6\cdot10^{-7}$ from the three pairs $\pm y,\pm2y,\pm3y$. The parent Fourier sign bounds cover all remaining vectors. Thus the Poisson summation denominator is at least $q(0)+\mu$. The resulting upper bound is $\max(T,U_n q(0)/(q(0)+\mu))$, which includes lattices of density at most $T$.

In rank 16, the parent's high-count bound is below $T=0.0172$, as are its frame bound and density thresholds. Above $T$ the count of short pairs therefore falls under the low-count proof. Its spatial bound is $p(0)+2(2^{16}-1)\kappa$, where $\kappa$ is the maximum of the archived class budgets. The replay reconstructs that expression and the parent's exact density bound, certifies $q>0$ by isolating its roots, and checks a contribution of at least $2.89\cdot10^{-12}$ from $\pm y$ over the full interval of possible squared lengths. Applying the same denominator estimate and taking the maximum with $T$ covers every count and every lattice. [Rank-16 check kit and scope](https://github.com/AstroNinja389/sphere-packings/releases/tag/dim16-check-kit-v6)

## Reciprocal planes in dimensions 19 and 23

The native proofs enclose a reduced dual plane by a rational box in three parameters. Write $a$ for the first squared length, $b$ for the squared length of the perpendicular projection of the second basis vector, and $\alpha\in[0,1/2]$ for its signed coefficient along the first vector after choosing signs. The squared lengths in Fourier coordinates are $a(i+\alpha j)^2+bj^2$ for integers $i,j$. The replay uses the same area, projection and reduction inequalities as the parents to discard impossible boxes. [Native plane domains and Fourier proofs](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/README.md)

Every remaining box is covered by lower bounds for the sum over a fixed finite set of opposite pairs. A mesh provides the first bound. In smaller boxes, the checker also encloses the endpoints and all real stationary points of the Gaussian polynomial on each possible length interval. The stationary points are roots of the exact polynomial $2q'-q$. This yields new plane contributions $\mu=7.25\cdot10^{-12}$ in rank 19 and $\mu=1.005\cdot10^{-12}$ in rank 23. If the old plane contribution was $\mu_0$, the new upper is $\max(T,U_n(q(0)+\mu_0)/(q(0)+\mu))$. The improvement is tiny because the first shells can lie close to a Fourier minimum in the reduced-plane domain.

## Joint spatial and radial-moment bound in dimension 36

Write $H(a)$ for the parent’s three compact radial features, $B$ for the constant harmonic block and $z$ for its origin vector. Replacing $z$ by $z+w$ changes the Schur cost by $2w^T B^{-1}z+w^T B^{-1}w$. It adds $2w^T H(a)$ to the spatial test at each nonzero lattice vector. The witness supplies $w$ and a decreasing list of rational upper bounds for that shifted test on successive squared-length intervals.

The verifier checks the interval heights with directed exponential arithmetic and uses the certified prefix counts to bound their sum. One opposite pair at squared length one is present in every normalized lattice. The verifier bounds its contribution at that exact point and reserves one place for it in each cumulative count. Both signs are included. The decrease in Schur cost exceeds the spatial allowance. Dividing that saving by the larger denominator $Q(0)$ gives a conservative density decrease on the parent’s high-density branch. Taking the maximum with its complementary threshold yields the universal bound. [Rational witness](joint-slack36.json)

## Replay

Use Python 3.11 or later with assertions enabled. From this directory, install the proof dependencies and supply the pinned `sphere-packings` checkout for 14D and the corresponding ZIP for each of the other dimensions. Keep the ZIPs intact for the archive hash checks. The rank-16 command uses the hashed upstream class-budget outputs. The upstream kit separately supplies a last-stage replay on Windows. The other commands run the corresponding parent replays and the underlying lower-rank exact-census inputs are documented by their proofs. [Pinned checkout](https://github.com/AstroNinja389/sphere-packings/tree/d7ffd0df47dfeb11814631858205f05ed73ec8ff) · [Native ZIP](https://github.com/AstroNinja389/sphere-packings/raw/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/native-proofs.zip) · [Rank-16 ZIP](https://github.com/AstroNinja389/sphere-packings/releases/tag/dim16-check-kit-v6) · [25–48D ZIP](https://github.com/AstroNinja389/sphere-packings/raw/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/n25-48-proofs-r9.zip)

```sh
python -m pip install python-flint==0.9.0 numpy scipy sympy mpmath
python -B verify.py 14 /path/to/sphere-packings
python -B verify.py 15 /path/to/native-proofs.zip
python -B verify.py 16 /path/to/n16-check-kit-v6.zip
python -B verify.py 19 /path/to/native-proofs.zip
python -B verify.py 23 /path/to/native-proofs.zip
python -B verify.py 36 /path/to/n25-48-proofs-r9.zip
```
