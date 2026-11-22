# Lattice sphere-packing density bounds

These certificates give strictly smaller upper bounds than the exact values in the benchmark for every dimension from 25 through 35 and for dimension 37. The bounds apply to all full-rank lattice packings in each listed dimension. The comparison is pinned to commit `d7ffd0df47dfeb11814631858205f05ed73ec8ff`. [Benchmark](https://github.com/AstroNinja389/sphere-packings/tree/d7ffd0df47dfeb11814631858205f05ed73ec8ff) · [Exact benchmark values](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/bounds-25-48.json)

| Dimension | Benchmark upper | New upper | Reduction |
|---:|---:|---:|---:|
| 25 | 0.001286558046689577 | 0.0012865566944499901305 | 0.000105% |
| 26 | 0.000869110293813701 | 0.0008690024388773178988 | 0.012410% |
| 27 | 0.000580593927514536 | 0.0005804875247073615033 | 0.018327% |
| 28 | 0.000391850205807269 | 0.0003917988242226201358 | 0.013113% |
| 29 | 0.000271857852999427 | 0.0002718573994663781847 | 0.000167% |
| 30 | 0.000189146062526746 | 0.0001891457703446318563 | 0.000154% |
| 31 | 0.000135429606979112 | 0.0001354139808594089906 | 0.011538% |
| 32 | 0.000095876900128018 | 0.0000958689943081682508 | 0.008246% |
| 33 | 0.000073452496873459 | 0.0000734453250082503981 | 0.009764% |
| 34 | 0.000052765505427207 | 0.0000527583769627631456 | 0.013510% |
| 35 | 0.000037469507068664 | 0.0000374650450434542721 | 0.011908% |
| 37 | 0.00001824556297305869842 | 0.0000182404738979538537 | 0.027892% |

The new decimals are rounded upward. The percentage is the reduction relative to the benchmark's exact rational upper bound.

## Replay

Use Python 3.11 or later with NumPy and python-flint. Download the [pinned proof archive](https://github.com/AstroNinja389/sphere-packings/raw/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/n25-48-proofs-r9.zip) and run the checker with assertions enabled. The checker accepts either the ZIP file or the directory extracted from it.

```sh
python -m pip install numpy python-flint==0.9.0
python -B verify.py /path/to/n25-48-proofs-r9.zip
python -B verify.py /path/to/n25-48-proofs-r9.zip 27 31 37
```

The replay is offline and requires no optimizer. It checks the archive hash and file manifest, replays the parent proofs, and checks each new bound against the exact benchmark value. The rows print as JSON. The benchmark archive supplies the original geometric and Fourier proof inputs. The new rational data are in `certificates.json`, `origin-certificates.json`, and `kernels/`.

## Angular counts in dimensions 25–35

Normalize the lattice to squared minimum one. For each source length window let $m$ count its opposite pairs of lattice vectors and let $m_j$ count the actual opposite pairs in the $j$-th difference-length bin. The parent proofs establish inequalities of the form

$$
b m^2 \le h m + \sum_j c_j m_j.
$$

Here $b>0$ follows from the degree-zero spherical coefficient, $h$ is the diagonal term, and $c_j\geq0$ pays for the allowed shorter signed differences. The replay checks those charges using the parent polynomials and the certified counts of decompositions per difference. The parents' matching inequalities with linear budgets are also available. [Parent angular and parity argument](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/notes/PROOF-parity-Fourier.md)

For any nonnegative integer $k$, the integer identity $(m-k)(m-k-1)\geq0$ gives a valid linear inequality

$$
((2k+1)b-h)m-\sum_j c_j m_j\le b k(k+1).
$$

Different tangents from several windows can be used at once because their statistics refer to the same lattice. For each nonnegative rational multiplier, add its constant contribution to the spatial allowance and subtract its weighted bin coefficients from the radial scores. Where grids cross, the smallest applicable coefficient gives a valid debit throughout the refined bin. Use the checked integer class profiles and the global prefix capacities to bound the shifted scores. The original spatial and Fourier signs are verified with the original Gaussian polynomials. The rank-31 proof also replays the benchmark's compact-kernel origin reallocation, which preserves the diagonal exactly.

## The compact kernel in dimension 37

Write the normalized kernel as a sum of angular harmonic matrices $A_j$ on the original radial features. Its diagonal depends on $\sum_j A_j$, and its origin contribution is $z^\mathsf{T}A_0^{-1}z/\alpha$. The new rational matrices are positive semidefinite and preserve $\sum_j A_j$ exactly. The checker proves the two-line and three-line inequalities over their full geometric domains by exact polynomial subdivision. The revised constant harmonic reduces the origin contribution, and the five-scale parent Fourier proof converts the saving into the density upper bound. [Benchmark compact-kernel proof guide](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/README.md)
