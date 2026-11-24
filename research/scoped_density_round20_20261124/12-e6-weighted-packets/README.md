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

## Joint medium/heavy bound

The supplied companion certificate `universal-dual.json` rules out a
set of 177 occupancy-at-least-four packets with 55 or more vertices
having no singleton-bad neighbors. In this graph, a packet of size at
least six is such a vertex: if two packet occupancies sum to more than
nine, the two phase sets realize every relative translation. In
particular this holds for a six-site packet against any packet counted
here.

Write $d$ for the singleton-bad degree of a uniformly chosen packet in
the 177-set. For integer $k$,
\[
\mathbb E[d^2]-2k\mathbb E[d]+k^2(1-55/177)\ge0
\]
follows from $(d-k)^2\ge0$ on the non-isolated vertices and $d=0$ on
at least 55 isolated vertices. The exact companion certificate uses
only the two degree cuts with $k=19,20$. No other bounds on the
singleton-bad degree are used.

Consequently, if there are at least 177 occupancy-at-least-four
packets, at most 54 of them have size at least six. Otherwise take
177 such packets retaining at least 55 of the size-at-least-six
packets. This implication permits arbitrary lighter packets elsewhere.

In a motif where every occupied packet has size 4, 5, or 6, let
$A,H,J$ be the number of occupied, size-at-least-five, and size-six
packets, respectively. The exact motif count is $4A+H+J$.
If $A\le176$, use $H,J\le108$ to obtain at most 920 sites.
If $177\le A\le189$, use $H\le108$ and $J\le54$ to obtain at
most 918 sites. Thus **any mixture of occupancies 4, 5 and 6 has at
most 920 motif sites**, below the 923 required. A prospective
improvement must contain a packet of occupancy 1–3 or 7–9. The
general equal-scale mask remains open.

Companion replay: `python -B verify_packets.py universal-dual.json`.
