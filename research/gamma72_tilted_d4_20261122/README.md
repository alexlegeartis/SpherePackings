# Tilted 77D and 78D lattice packings

Both lattices have squared minimum 8. Tilting the new generators into the old
D4 directions reduces the Gram determinant without reducing the minimum.

| Dimension | Gram determinant, old → new | Occupied density at least | Gain over the preceding repository bound |
|---:|---:|---:|---:|
| 77 | 8 → 7 | 6.23325752266126e-16 | 6.904496% |
| 78 | 12 → 10 | 2.08655477067461e-16 | 9.544511% |

These improve the preceding [lattice and unrestricted lower bounds](../gamma72_depth6_20261122/README.md).
The comparisons to the named Bachoc–Nebe sections are separately tabulated in [results.json](results.json).
These are explicit examples, not optimality claims or an exhaustive literature priority claim.

## Why the tilt works

Keep the verified 76D extension of [Nebe’s Γ72](https://arxiv.org/abs/1008.2862) by standard
D4. Each new direction projects to e₁/2 on D4. Give the residual one-dimensional
Gram entry 7/4 in 77D and the residual Gram
`[[7/4, -3/4], [-3/4, 7/4]]` in 78D. The determinants are respectively
4·(7/4) = 7 and 4·(5/2) = 10. Orthogonal A1/A2 attachments to D4 had larger determinants.

The horizontal Γ72 shifts are specified in [phase.json](phase.json) and
[phases.json](phases.json). Each isolated shift has squared clearance at least 6.
For every nonzero new parity, including the sum of the two shifts in 78D,
all 16 old D4 parity translates have two independent A labels and no
pair of norm-four Leech blocks with the required total B label.
This excludes horizontal squared norm below 5. A single-layer step has
vertical squared norm at least 3 except at the unshifted D4 point, where its
vertical squared norm is 2 and the stronger depth 6 applies.
In 78D a step with both coefficients 1 has residual squared height 2,
and its projection e₁ has squared distance 1 to D4 because its coordinate
sum is odd. The verifier covers all other short coefficient pairs.
The embedded Γ72 attains squared norm 8.

## Replay

The exact Grams are `gram_numerator/4` in [gram77.json](gram77.json)
and [gram78.json](gram78.json). Run with Python, NumPy and SymPy:

```sh
OPENBLAS_NUM_THREADS=1 python verify77.py
OPENBLAS_NUM_THREADS=1 python verify78.py
OPENBLAS_NUM_THREADS=1 python verify_short_coefficients.py
```

The first two verifiers replay the pinned predecessor, independently check
the new horizontal syndromes and rebuild the literal Grams. The third
enumerates the complete coefficient box given by the inverse vertical Gram,
checks all 105 short tilted parity/norm cases in 78D (including the 77D
sublattice), and verifies the published density floors and comparisons.
