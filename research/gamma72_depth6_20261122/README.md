# Six explicit Γ72 lattice density improvements (73–78D)

| Dimension | Occupied density at least | Gram determinant | Improvement over previous repository lattice |
|---:|---:|---:|---:|
| 73 | 4.26462994945094e-14 | 2 | 11.8034% |
| 74 | 1.43007109517410e-14 | 3 | 25% |
| 75 | 5.05260043732660e-15 | 4 | 39.7542% |
| 76 | 2.04778582360778e-15 | 4 | 56.25% |
| 77 | 5.83067851333243e-16 | 8 | 46.1585% |
| 78 | 1.90475519228084e-16 | 12 | 36.71875% |

All are lattices of minimum squared norm 8. The 78D entry also improves
the preceding repository unrestricted packing 1.41519803035180e-16
(approximately 34.59% gain). The lattice comparison determinants are
5/2, 75/16, 125/16, 625/64, 4375/256, and 91875/4096, respectively,
from the preceding [Γ72 construction](../gamma72_leech_layers_20261122/README.md).
This is a comparison to recorded constructions, not an assertion of world priority or unrestricted optimality.

## Why the half cosets are deep

Use Nebe’s [Γ72 theorem](https://arxiv.org/abs/1008.2862) and
[catalogue Gram](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/neb72.html).
In normalized block coordinates, Γ72 is a sublattice of three Leech lattice
blocks with one-half of their direct sum metric. Its image consists exactly
of triples whose mod-2 common components on one polarized half A agree
and whose complementary B components sum to zero. These are 36
independent binary constraints; the exact Gram identity verifies the index.
An integral block-coordinate shift is therefore a half coset of Γ72.
It is identified by two binary A syndromes (relative to block 3)
and the total binary B syndrome.

A shift with independent A syndromes a,b cannot be represented in a
single occupied Leech block. To rule out total squared Leech norm below
12, it remains to exclude two blocks of norms 4+4 or 4+6. The possible
pairs of nonzero A labels are exactly the three distinct pairs in
{a,b,a+b}. In the exhibited A triangle, all norm-four B-lists
annihilate the whole triangle in the mod-2 Leech bilinear form; the
chosen total B syndrome does as well. A norm-six B label would therefore
annihilate its own A label, impossible: in the Leech quadratic form
modulo 2, a norm-six class has quadratic value 1, while both
polarization halves are totally singular. An exact enumeration of the
196,560 norm-four vectors then rules out all remaining 4+4 B-sums.
This gives depth at least 6 in Γ72's halved metric.

A nonzero B-only syndrome has depth at least 4: a single block of norm 4
would require a norm-four vector with A label zero, which is absent in
the enumerated shell, and a single block of norm 6 violates singularity.
More generally, two independent A syndromes require two occupied blocks,
and if the A syndromes are nonzero and dependent the verifier excludes
norm 4 using the exact list and norm 6 using the orthogonality of the
required B label with the nonzero A label. Every nonzero syndrome
automatically has squared depth at least 2.

## Vertical lattice and replay

The 73D vertical lattice is A1. For 74–76D the lifts use leading A2, A3,
D4 Cartan blocks. The root parities of D4 map onto a two-dimensional
space of A-syndrome pairs in which all three nonzero matrices are
invertible on the selected A triangle. The three nonzero non-root
parities instead have B-only colors. For 77–78D the vertical lattice
is D4+A1 or D4+A2. The last one or two directions use the opposite
orbit of invertible A-syndrome pairs. Every required short parity,
including cross-block norm-four cases, is checked separately.

If G is the original Gram, W is the integer numerator matrix of the
half shifts, and V is the vertical Gram, the exact full Gram is

```
 [ G          G W^T/2           ]
 [ W G/2      W G W^T/4 + V      ] .
```

Its determinant is det(V) by the Schur complement. For vertical
coefficients z with norm squared below 8, the componentwise exact
bounds \(z_i^2\leq 8(V^{-1})_{ii}\) make the minimum check finite;
norm 2, 4, and 6 coefficients use the corresponding depths above.
The embedded Γ72 attains norm 8.

The [explicit numerators and vertical Grams](numerators.json) specify all six lattices.
For each one the verifier reconstructs the full Gram and checks its determinant.
[Exact density squares and comparisons](results.json) distinguish the preceding
repository lattices, the preceding unrestricted constructions, and the named
external comparisons. In 78D the named external comparison is the L80
section of determinant 48 with squared minimum 8.

To replay the checks, run
`OPENBLAS_NUM_THREADS=1 python verify.py` with NumPy and SymPy.
The verifier checks the occupied-density floors with rational Machin bounds.
