# High-clearance periodic five-torsion motifs in the Kappa11 torus

No new density lower bound. Consider any finite motif \(C\) in \((1/5)K/K\), with period \(K=\mathrm{Kappa11}\) of Gram determinant 972 as in the adjoining lattice-only proof. Do not assume that \(C\) is additive, or scalar invariant. If the squared minimum of the periodic set is \(d/25\) for an even \(d\ge38\), its center density is at most the known general 11D lower bound \(9/256\). The latter is the center-density normalization of the dimension-11 comparison in Cohn's [density table](https://cohn.mit.edu/sphere-packing/) and the [three-point bounds paper, Table 1.1](https://arxiv.org/abs/2206.15373).

The exact comparison used is
\[
 65536 |C|^2 d^{11}\ \le\ 972\cdot81\cdot100^{11}.
\]
In contrast to additive subspaces, a motif's differences can be compatible even if their doubles are not. Thus use the true least lift norm \(m(x)\), without scalar saturation, and restrict the four-isometry action to add only inversion. There are 1416 point/character orbits each. A sorted exact LDL search recomputes the least lift norm of every point representative. For a character orbit \(O\), let \(c=|O|\) and \(h_j\) be the number of characters taking pairing value \(j\bmod5\) at a representative point. Inversion gives \(h_j=h_{5-j}\), so four times its character sum is
\[
 (5h_0-c)+2(h_1-h_2)\sqrt5.
\]
At each tested threshold, use the supplied nonnegative integer weights on nontrivial character orbits and sum over all ordered motif differences, as in the lattice-only proof. The replay strictly bounds the function above by a negative constant on all nonzero allowed differences through rational bounds on \(\sqrt5\). It recomputes the pairing statistics from orbit inputs; no optimizer or saved Fourier matrix is used by the checker.

| Raw squared minimum d | Upper bound on size | Smallest size to beat 9/256 |
|---:|---:|---:|
|38|196|225|
|40|158|170|
|42|104|130|
|44|81|101|
|46|52|79|
|48|41|63|
|50|23|50|
|52|14|40|
|54|13|33|
|56–64|11|at least 13|
|66–100|1|at least 2|

The last line follows because the maximum least norm among nonzero classes is exactly 64. The supplied certificates at \(d=42,56\) reproduce earlier receipts, while the other listed thresholds below 56 complete the previously open high-clearance range. None of the thresholds \(d<38\) is ruled out by this note; in particular raw 34 and 36 remain possible in the pair-distribution relaxation used in this scan. No conclusions about arbitrary real phases, other periods, or global 11D sphere packing follow.

## Replay

    python check_inputs.py
    g++ -O3 -std=c++17 short.cpp -o short
    ./short 4 0 0 < K11.ldl
    g++ -O3 -std=c++17 fullhigh.cpp -o fullhigh
    ./fullhigh
