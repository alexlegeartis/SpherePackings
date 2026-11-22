# Lattice packing upper bounds in dimensions 38–48

The certificates give strictly lower density upper bounds in all eleven dimensions than the exact values in the pinned AstroNinja389 archive. They also give bounds 5.13% to 23.66% below Nina Zubrilina’s separate Gaussian-certificate results. Dimensions 41 and 42 use the values in her companion manifest. [AstroNinja389 benchmark](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/bounds-25-48.json) · [Zubrilina, Table 1](https://arxiv.org/html/2610.35830v1#S4.T1) · [Companion manifest](https://github.com/AstroNinja389/lattice-packing-bound/blob/340a0e93bb95942c2b1b80092be43149504970eb/manifest.json)

| Dimension | New lattice density upper bound | Reduction vs. AstroNinja389 | Reduction vs. Zubrilina |
|---:|---:|---:|---:|
| 38 | 0.0000123421888282022774971681 | 0.005801033% | 19.1536% |
| 39 | 0.0000084768304132050218168335 | 8.055632e-05% | 22.5767% |
| 40 | 0.0000059759954747923372601960 | 2.515234e-08% | 23.6603% |
| 41 | 0.0000043915394964222208569670 | 9.804743e-08% | 20.5974% |
| 42 | 0.0000034723546103562322281453 | 1.455467e-08% | 10.4639% |
| 43 | 0.0000025359881069348761410667 | 9.710723e-08% | 6.1452% |
| 44 | 0.0000017667499405165497949257 | 2.000000e-13% | 5.1311% |
| 45 | 0.0000011727225939886516583185 | 2.405196e-06% | 8.4057% |
| 46 | 0.0000007717674119868513114249 | 1.631165e-06% | 12.3665% |
| 47 | 0.0000005127256280869701134112 | 1.020594e-06% | 15.5487% |
| 48 | 0.0000003453056409988827388707 | 4.126924e-07% | 17.7742% |

Density decimals are rounded upward. Percentages are relative decreases computed from the exact fractions in [certificates.json](certificates.json). The largest new gain over AstroNinja389 is about 0.0058%. The [research notes](RESEARCH.md) develop approaches for a 1% reduction against that benchmark.

The benchmark bound table and archive are unchanged through the repository head checked on November 22, 2026. [Checked head](https://github.com/AstroNinja389/sphere-packings/tree/7b63ee2094b7489aea3da83349c547c8da4cb388)

## Offline replay

Use Python 3.11 or later with NumPy and python-flint, and download the [pinned proof archive](https://github.com/AstroNinja389/sphere-packings/raw/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/n25-48-proofs-r9.zip). The verifier accepts the ZIP or its extracted directory. Every run checks the parent manifest and replays the applicable parent proofs.

```sh
python -m pip install numpy python-flint==0.9.0
python -B verify.py /path/to/n25-48-proofs-r9.zip
python -B verify.py /path/to/n25-48-proofs-r9.zip 38 42 44 48
```

Append `--nina-manifest /path/to/manifest.json` to compare the Zubrilina column against the pinned companion manifest. Its eleven certificates also replay with the companion repository’s `verify/run.py --verifier bernstein` and the selected files `certificates/n038.json` through `certificates/n048.json`.

The new sign and matrix checks use exact rational arithmetic. Ball arithmetic encloses the volume factors and exponential terms. The verifier runs with assertions enabled and requires no optimizer.

## Angular energy in dimensions 38–43 and 45–48

As in the parent proofs, work with a densest lattice and normalize its squared minimum to one. In a compact kernel let $H(a)$ be the vector of radial features at squared length $a$, $B_j$ the harmonic blocks after restriction to dimension $n$, and $z$ the origin vector in the same normalization. The parent proof uses $s_0=z^{\mathsf T}B_0^{-1}z$ as the allowance for the origin. Its off-diagonal and radial inequalities continue to hold for the original kernel.

The witness selects a radial direction $v$ with $w(a)=v^{\mathsf T}H(a)\ge0$ and an angular polynomial $F(t)=\sum_j c_j G_{2j}(t)$ with $c_j\ge0$. Exact principal minors give $B_j-c_jvv^{\mathsf T}\succeq0$. If $a,b\in[1,C]$ for two distinct support lines, the squared minimum gives $|t|\le(a+b-1)/(2\sqrt{ab})\le 1-1/(2C)$ since $C\le2$. The verifier checks $F$ on that interval or proves it is a positive multiple of a square. It also proves the pointwise radial inequality $w(a)^2\ge\ell^{\mathsf T}H(a)$.

For $u=2\sum_{\text{lines}}H(a)$, the angular subkernel contributes at least $4F(1)\sum_{\text{lines}} w(a)^2\ge2F(1)\ell^{\mathsf T}u$. The residual kernel contributes at least $u^{\mathsf T}(B_0-c_0vv^{\mathsf T})u$ by harmonic positivity. With the original linear term $2z^{\mathsf T}u$, completing the square gives the new allowance

$$
s_1=(z+F(1)\ell)^{\mathsf T}(B_0-c_0vv^{\mathsf T})^{-1}(z+F(1)\ell).
$$

The difference $s_0-s_1$ is positive in each recorded row. Subtract its downward-rounded density contribution from the exact parent upper bound. In ranks with a conditional parent bound, keep the same density threshold and denominator contribution.

## Positive dual mass in dimension 44

Write $U$ for the exact parent density bound, $U_{n-1}$ for its predecessor, $V_k$ for the unit ball volume, and $S=4\pi^2/\tau$ for the Fourier squared-length scale. Suppose the densest lattice has density above $T=0.999U$. A section perpendicular to a primitive shortest dual vector $y$ has determinant $\det L\,\|y\|$. The section bound gives

$$
S\|y\|^2\ge S\left(\frac{2V_{n-1}T}{V_nU_{n-1}}\right)^2.
$$

Hermite's inequality on the scaled dual gives the upper end of the interval. On the certified interval $[a,b]$, rational Bernstein checks on intervals of width at most one give $Q(x)e^{-x/2}\ge\eta>0$. The original dual polynomial is globally nonnegative. The pair increases the Poisson denominator from $D$ to at least $D+2\eta$, which gives the bound $\max(T,UD/(D+2\eta))$.
