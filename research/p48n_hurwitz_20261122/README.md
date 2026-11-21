# Explicit quaternionic lattice packings in dimensions 51 and 52

| Dimension | Occupied density at least | Increase over previous repository bound | Gram determinant |
|---:|---:|---:|---:|
| 51 | **9.34489307985541e-10** | **17.2603%** | 4 |
| 52 | **3.95931817049999e-10** | **37.5%** | 4 |

These lattices have squared minimum 6. Their center-density squares are exactly \((3/2)^n/4\). The comparisons use the exact [preceding P48q extension certificates](../adaptive_layers_20261115/P48_PROOF.md). The repository's [dated literature review](../covering_and_deep_cosets_20261115/PRIORITY.md) documents the earlier external comparators.

## Construction and proof

Take the even unimodular 48D lattice \(P_{48n}\) of squared minimum 6 from Nebe's [catalogue entry CQ48a](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ48a.html). Its published minimum is the parent theorem. The [original catalogue snapshot](../density_breakthroughs_20261113/seed_deformations/sources/CQ48a.html) supplies its exact Gram \(G\).

The certificate gives integral quaternion operators \(I,U,X,Y\), where \(X^2=Y^2=-I\), \(XY=-YX\) and \(2U=I+X+Y+XY\). It also gives a binary seed \(w\) in a specified reduced basis. A complete enumeration excludes all vectors of squared norm at most 12 in \(w+2P_{48n}\). Their norms are multiples of four, so the corresponding horizontal half-coset has squared depth at least 4.

Equip the vertical space with Gram
\[
M=
\begin{pmatrix}
2&1&0&0\\
1&2&1&1\\
0&1&2&0\\
0&1&0&2
\end{pmatrix}.
\]
For the 51D packing use its leading rank-three block. For each vertical coefficient \(z\), take the horizontal shift specified by half the numerator \(w(z_1 I+z_2 U+z_3 X+z_4 Y)\), expressed in the catalogue basis. The verifier enumerates the 24 vertical roots and checks that every corresponding operator is an integral isometry of \(G\). All twelve root parities thus inherit depth 4 from the single seed. In the remaining three nonzero parities, the horizontal class is nonzero and has numerator norms divisible by four. The parent minimum excludes norms below 6, so these three horizontal depths are at least 2. Their vertical squared norms are at least 4 since the even vertical lattice has no roots in those parities. All-even nonzero vertical coefficients already give squared norm at least 8. These cases give squared minimum 6, which occurs in the embedded parent.

The exact numerator rows \(B\) are in [results.json](results.json). In a lattice basis with the 48 parent coordinates first, the Gram is
\[
\begin{pmatrix}
G & GB^T/2 \\
BG/2 & BGB^T/4+M
\end{pmatrix}.
\]
Its determinant is 4 in both dimensions. The verifier also checks the distinctness of the 16 horizontal parity classes, the unimodular basis change, the exact determinants and the outward density rounding.

## Replay

Install python-flint and SymPy. Run `python research/p48n_hurwitz_20261122/verify.py` to verify the exact setup, archived enumeration reports and comparisons. Add `--replay` to compile and rerun the full seed exclusion and the short-vector control. The full replay uses GNU C++ with IEEE-754 binary64 and the existing outward-rounded [interval enumerator](../covering_and_deep_cosets_20261115/explicit_p48/interval_enum.cpp). Every leaf norm is checked with integers.
