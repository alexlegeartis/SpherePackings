# Occupancy-four packet obstruction for 12D ternary E6×E6 at equal scales

We consider only the equal-scale ternary motif mask on `V×V`, where
`V=E6/3E6`, with separation requiring the sum of the two least squared
integer-lift norms to be at least 8. Let `R` be the ternary radical of
order three, and group motif points into packets on cosets of `R×R`.
The claim concerns **all packets having at least four sites**, allowing
arbitrary remaining packets of size at most three.

The rational three-point certificate in `packet-dual.json` excludes
190 such packets and hence every larger collection by taking a subset.
Thus there are at most 189 packets with at least four sites. In particular, it excludes an
improvement over Coxeter–Todd by an equal-scale motif in which every
occupied packet has exactly four sites: that would require at least
231 packets (923 sites). This is not a global density bound. See the
end for mixed-weight consequences.

## Affine forbidden offsets and triangle holonomy

The rank-five quotient `W=V/R` has shell types 0 (zero), 1 and 2 (the
two nonzero quadratic values), and 3 (nonzero isotropic). The three
radical lifts of a displacement of each respective type have squared
depths listed as multisets
`(0,12,12), (2,8,8), (4,4,4), (6,6,6)`.
This is replayed by integer norm enumeration, using the Cauchy bound
on each coordinate derived from the exact inverse Gram matrix.

For distinct packet labels there are three cases:

- **Good** types `(1,3),(2,2),(2,3),(3,1),(3,2),(3,3)`: no forbidden
  radical difference in the nine-site affine plane.
- **Singleton-bad** types `(0,1),(1,0),(1,1)`: exactly one forbidden
  radical phase difference.
- The other six off-diagonal types: forbidden differences form an
  affine line in the radical plane.

Every four-point subset of the radical plane projects onto at least
two points modulo any line, since an affine line has only three points.
Two such projected sets realize every relative difference in the
three-point quotient. Hence two packets of size at least four cannot
have an affine-line forbidden condition.

A **triangle with all three edges singleton-bad** is also impossible
for three such packets. The unique forbidden offset in each block is
the minimal-depth lift of the corresponding displacement; it is zero
for type 0 and the unique root residue for type 1. These lifts are
additively coherent on a triangle involving only types 0 and 1.
For the nonparallel case, two norm-2 root lifts with type-1 difference
have inner product 1 modulo 3. Their inner product is in
`{-1,0,1}`, so it is 1 and their difference is itself a norm-2 root.
Zero and proportional cases hold directly. Thus the offsets around the
triangle sum to zero and the three packet phase sets can be shifted into
a common radical plane so that the forbidden conditions require pairwise
disjointness there. Their total occupancy, at least twelve, exceeds
nine. The independent replay verifies this coherence on every pair of
rank-five residues concerned, not just on a sample.

The same affine-line exclusion applies to any two noncollinear
three-point packets, and the all-singleton-bad triangle exclusion
holds for any three such noncollinear shapes if their total size
exceeds nine (thus it does **not** apply to three size-three packets).

## Exact separating certificate

The quotient translation group stabilizer orbitals on `W×W` use
`(type(x),type(y),type(x−y),x=−y)` in each coordinate. The Gram matrix
and linear dependency information determine these orbits by Witt's
extension lemma. All vertex permutation actions on these orbitals
are checked exactly. Taking products of the 40 factor orbitals and
symmetrizing by vertex permutations and the factor swap gives 89
nonnegative variables after deleting forbidden pairs and triangles.

If $C$ were a quotient code of size $N$, define $t(x,y)$ by
averaging the rooted product of three site indicators over all $N$
possible translating roots in $C$ and over the stabilizer.
Then $t(x,y)$ and its complement $a(x-y)-t(x,y)$ are positive
semidefinite as sums of outer products with the translating root
inside or outside $C$, respectively; $a(x)=t(0,x)$ is the averaged
pair autocorrelation. The root diagonal is $t(0,0)=1$ and the sum
over $y$ of $t(x,y)$ is $N t(0,x)$. These identities use exact integer
orbital multiplicities.

Each Rayleigh vector is an integer linear combination of products of
ordinary 243-entry integer vectors; it is not assumed to lie in a
numerically computed invariant block. The provided nonnegative
rational Rayleigh multipliers and unrestricted equality multipliers
give a separating linear functional with right-hand side 1 and strictly
negative coefficients of all 89 nonnegative variables.

Replay: `python -B verify_phases.py && python -B verify_packets.py`
(Python and NumPy only; exact integer and rational arithmetic).

## Combining the two occupancy thresholds

Let `A` be the number of packets of size at least four, and `H`
the number of size at least five. This certificate bounds `A` by
`N−1`. The independently verified heavy-packet bound gives
`H≤108`. Therefore if all occupied packet sizes lie in `[4,s]`,
the total motif count is at most `4(N−1)+(s−4)108` for `s≥4`.
More generally the same inequality bounds just the contribution of the
size-at-least-four packets in an otherwise arbitrary motif. It does
not bound the contribution of lighter packets.


## Round 21: all mixtures of occupancies four through seven

**Theorem.** If every occupied packet has size in `{4,5,6,7}`, the
total number of motif positions is at most 922, below the 923 required
to beat Coxeter–Todd. No arbitrary mixed motif bound is claimed.

A universal-good label is one adjacent in the good mask to every
other selected label. Occupancy at least six forces this property
among occupancy-at-least-four packets, since the two phase sets
have sizes summing to more than nine and therefore realize every
relative radical-plane difference.

For a hypothetical code of size N with at least h universal-good
labels, let d be the degree of a uniformly chosen label in the
singleton-bad graph. The function (d-k)^2 is nonnegative on the
other labels, and d=0 on at least h labels. Therefore
E[d^2] - 2k E[d] + k^2 (1 - h/N) >= 0.
The expected degree and its square are exact linear forms in the
rooted pair and triple variables, respectively. The five certificates
below use only these degree inequalities with k between 19 and 21,
plus ordinary three-point PSD constraints and row-sum identities.
Every numerical ingredient was replaced with exact integer Rayleigh
vectors and rational multipliers.

The certificates rule out the following (N,h) pairs. Each applies to
all codes of size at least N and with at least h universal-good
vertices, by taking a subset retaining h of them.

| Actual occupied packet count A | Excluded (N,h) | Upper bound on size-at-least-six count J |
|---|---:|---:|
| 150–159 | (150,88) | 87 |
| 160–169 | (160,68) | 67 |
| 170–176 | (170,54) | 53 |
| 177–183 | (177,40) | 39 |
| 184–189 | (184,30) | 29 |

Let H be the number of packets of size at least five. The included
good-mask certificate proves H <= 108 (two heavy packets are good
to one another). Also A <= 189 by the occupancy-four structural
certificate. Since sizes are at most seven, the literal motif
count is at most 4A+H+2J. For A <=149, use H,J <=108 to get at most
920 positions. In the five other bands, this gives bounds 918,
918, 918, 918, and 922, respectively, by the table. This proves
the claim.

Replay with Python and NumPy:
```sh
python -B verify_phases.py
python -B verify_heavy.py
python -B verify_packets.py
for p in vcert*.json; do python -B verify_packets.py "$p" || exit 1; done
```

## All mixtures of occupancies four through eight

**Stronger theorem.** Any motif in which every occupied packet has
size in `{4,5,6,7,8}) has at most **922** positions. This rules out
a Coxeter–Todd improvement within every mixture of those sizes.

The marked three-point certificates explicitly distinguish a color U
of h universal-good vertices from the remaining color L. If there
are more than h universal vertices, choose any h for the mark.
Edges involving a marked vertex must be good whenever endpoints
differ, and a triangle with all three edges singleton-bad remains
forbidden regardless of colors. Coincident vertices have the same
color. There are 372 nonnegative orbit variables under the same
stabilizer, factor swap, and simultaneous triangle/vertex-color
permutations.

For each root color r, let the rooted three-point matrix have
entries indexed by an offset and an endpoint color on each side:
`T_r[(x,c),(y,d)]` is the stabilizer-average of
`(1/N) sum_{z in C_r} 1_{z+x in C_c} 1_{z+y in C_d}`.
The matrices `T_L,T_U` and the complement with the root outside both
colors are positive semidefinite, as sums of outer products. The
unconditioned autocorrelation entry for colors c,d at difference
y-x is the corresponding rooted entry with the first offset zero
and root color c. The identity `sum_y T_r[(x,c),(y,d)] =
|C_d| T_r[(0,r),(x,c)]` holds pointwise; the verifier sums it by
exact displacement shells. The two origin entries are fixed to
`|C_L|/N` and `|C_U|/N`. The rational witness for each target
combines these exact equalities with 135–163 positive integer-vector
Rayleigh inequalities. All 372 resulting coefficients of nonnegative
variables are strictly negative, while the right-hand side is
strictly positive.

| Actual occupied packet count A | Excluded marked (N,h) | Upper bound on size-at-least-six count J | Bound on total at sizes 4–8 |
|---|---:|---:|---:|
| 123–139 | (123,87) | 86 | 922 |
| 140–149 | (140,73) | 72 | 920 |
| 150–159 | (150,60) | 59 | 921 |
| 160–169 | (160,47) | 46 | 922 |
| 170–179 | (170,33) | 32 | 920 |
| 180–189 | (180,20) | 19 | 921 |

Use the same A,H,J counts as above. The total is bounded by
`4A+H+3J`, where `H<=108` and `A<=189`. If `A<=122`,
then `4A+4*108<=920`. Otherwise, take an N-subset retaining h
universal-good vertices if present and apply the marked target in
the table. Thus a record in the equal-scale mask now requires an
occupancy-1–3 packet or a full nine-point packet.

For exact marked replay (Python and NumPy only):
```sh
for p in mcert*.json; do python -B verify_marked.py "$p" || exit 1; done
```

## Full packets allowed: exclusion at high packet count

The same marked method also rules out a record-size contribution
from packets of sizes **4–9** if their number A is at least 118.
This conclusion still allows an arbitrary population of lighter
packets elsewhere, but does not bound their contribution.

For this submotif the total is at most `4A+H+4J`. Thus a contribution
of at least 923 requires `A+J>=204`, via `H<=108`. The following
exact marked targets exclude all such cases with A at least 118
(each marked exclusion `(N,h)` excludes `A>=N, J>=h)):

| A | (N,h) |
|---|---:|
| 118–120 | (118,84) |
| 121 | (121,83) |
| 122–124 | (122,80) |
| 125–128 | (125,76) |
| 129–132 | (129,72) |
| 133–137 | (133,67) |
| 138–140 | (138,64) |
| 141–147 | (141,56) |
| 148–154 | (148,49) |
| 155–161 | (155,43) |
| 162–168 | (162,36) |
| 169 | (169,30) |
| 170–176 | (170,28) |
| 177–182 | (177,22) |
| 183–189 | (183,15) |

Thus any record-size submotif made only of packets above size three
must have between 103 and 117 occupied packets, including at least
87 with at least six positions, and at least 23 full nine-site
packets. For the last claim, replace the occupancy of every full packet
from nine to eight; with F full packets the original count is at most
`4A+4H+F <= 4*117+4*108+F = 900+F`.
The lower limit 103 is the raw nine sites per packet.
This does not rule out the low-count region.

## Fully radical-invariant additive motifs and affine rank-three jams

The good-mask quotient has no linear code of dimension four. Its
projections would both be injective, so the induced quadratic forms on
a four-dimensional code subspace each have rank at least three
(a four-dimensional subspace of nondegenerate rank five has radical
dimension at most one). The exact program `verify_additive.py` tests
the four rank/discriminant possibilities for the first form. Where it
takes value 1, the second form must vanish; where it takes value 2,
the second form must avoid value 1. In all four cases these tests
force the second form to be identically zero, a contradiction.
So an additive full-packet motif has at most 27 packet labels and
243 sites. An explicit diagonal rank-three code attains this
within the restricted family.

Moreover, **every affine rank-three good code is maximal**, even
against arbitrary nonlinear good additions. The exact
`verify_linear_jam.py` first normalizes one induced rank-three-domain
quadratic form to a diagonal form, enumerates all 729 homogeneous
quadratic forms for the second projection, and finds 30 ordered
form cases with injective realizations. Witt extension on each
five-dimensional factor ensures that the pair of normal forms and
one injective realization per factor cover all embeddings including
degenerate subspaces. The code tests every possible new label for
every resulting embedding; in all cases the common good neighborhood
is empty. Affine translates reduce to the linear case.

In particular, a rank-three affine subspace of full packets blocks
every additional occupied packet, even a light one: nine plus any
positive occupancy realizes all radical offsets, so every new label
would have to be in the nonexistent common good neighborhood.
```sh
python -B verify_additive.py
python -B verify_linear_jam.py
```

Extract the adjacent replay.tar.gz, then run the documented
commands from the extracted release/ folder. The archive is
the full exact source and witness data, not numerical solver output.
