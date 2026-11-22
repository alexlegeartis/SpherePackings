# MW92 projected-simplex packings in 86D and 87D

| Dimension | Translates | Occupied density at least | Over preceding unrestricted bound |
|---:|---:|---:|---:|
| 86 | 7 | 3.2157847626326119e-19 | +13.0418% |
| 87 | 6 | 2.6231014000959713e-19 | +1.6557% |

These are explicit periodic packings, not lattice-only bounds or optimality
claims. Relative gains are rounded down and use the preceding [86D/87D
correlated MW92 construction](../packing_release8590_20261117/).
The external comparison in both dimensions is the
[derived 1982 Bos–Conway–Sloane section](../priority_resolution_20261115/CLASSICAL_SECTIONS.md).

## Proof

Let `G` be the Gram of Shioda's MW92 lattice in the
[pinned parent certificate](../mw92_sections_20261115/). It has minimum
squared norm 16 and determinant 47^14, as in the
[Schütt–Shioda treatment, Chapter 13](https://doi.org/10.1007/978-981-32-9301-4_13).
For each result, the corresponding witness stores a rank-k integer matrix
`B` and a list of distinct labels in Z^k, including zero, where k=92−n.

The exact column HNF of `B` is the identity, so the pairing map from the
parent coefficient lattice Z^92 to Z^k is onto. Let `H=B G^−1 B^T`. The
normal component of a lift of label u is `G^−1 B^T H^−1 u`. In
particular the squared normal length lost when projecting a difference
of labels u,v is exactly `(u−v)^T H^−1 (u−v)`.

For each label, project all parent vectors with that label onto the
orthogonal complement of the normal space. This gives one translate
of the rank-n kernel lattice. In a translate, differences keep the
parent minimum. A pair from different labels comes from two distinct
parent vectors, so its projected squared separation is at least 16−β,
where β is the maximum of the normal loss over all label pairs.
It is positive in both dimensions. This also proves distinctness for
all translates and checks all cross-coset translates, not just differences
between chosen representatives.

The normal projections for all integer labels form a full lattice
with Gram `H^−1`. Coarea gives kernel squared covolume 47^14 det(H).
Thus the occupied density is at least the unit-ball volume times the
number of labels times `((16−β)/4)^(n/2) / sqrt(47^14 det(H))`.
The exact values of β, period determinant and squared center density
are in [verified_results.json](verified_results.json).

The search used dual normals of squared norm 182/47 and targeted
near-simplex correlations between them. Verification uses only the
literal normals and labels and does not rely on the search.

## Replay

From the repository root, use Python with NumPy and SymPy:

```
OPENBLAS_NUM_THREADS=1 python research/mw92_small_simplexes_20261122/verify.py
OPENBLAS_NUM_THREADS=1 python research/mw92_small_simplexes_20261122/audit.py
```

The second script reconstructs the inverse in exact arithmetic and
rechecks the coarea and projection calculations from the lattice lifts.
