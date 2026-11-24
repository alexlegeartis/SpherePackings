# Exact rigidity of the BN80 + D4 84D packing

This is a scoped result about the already verified 84D incumbent, not a new
density lower bound or an assertion of global optimality. The incumbent has
squared minimum 8 and covolume 2. Write H for the minimum-eight integral even
unimodular rank-80 Bachoc–Nebe lattice and D for the vertical D4 with its simple
root Gram V. The packing projects onto D, with full translates of H above each
v in D. Its horizontal phase map F: D -> (1/2)H/H depends on the four root
coefficient parities. The input is the verified BN80 / D4 root-basis chart.

## Literal and finite data

The independent exact replay verifies all 24 vertical roots, that their outer
products span Sym4, and that their outer product sum in Euclidean coordinates
is 12 I4. It checks every one of the 12 root colours has a norm-6 vector in its
horizontal phase. For each of the four simple root colours (masks 1,2,4,8), it
checks 80 such vectors with rank 80 modulo 101 and hence over R. These are
actual contact vectors because the packed minimum is 8 and the associated
vertical root has norm 2. Their negatives are in the same phase because it
has order two. It also verifies a norm-6 witness for each of the three
non-root, nonzero colours. No nearest-vector search is trusted by this replay.

An immediate global scoped consequence: if a lattice extension has horizontal
intersection H and all its horizontal phases are restricted to these 16 fixed
labels, then its density is at most that of the incumbent. Indeed its squared
minimum mu is at most 8. When its full lattice projects to A, the finite
phase image (order at most 16) ensures A is discrete. Its kernel A0 in the
vertical space has squared minimum at least mu and index at most 16 in A.
The exact gamma_4=√2 gives covol(A)≥mu^2/(2*16), hence centre density at most
2(mu/4)^40≤2^41. The same proof applies to arbitrary phases with image size
at most 16. At mu=8 our exact 15 depth witnesses additionally force vertical
separation at least sqrt(2) for every full-fibre pairing using these labels
(even nonadditive motifs). This last projection does not assert unrestricted
D4 optimality among arbitrary 4D packings.

## Strict local optimum under arbitrary lattice strains

The norm-8 shell of H is a spherical 4-design. Indeed its degree-2 and
degree-4 harmonic theta series have respective weights 42 and 44 and order
at least 4 at infinity, impossible for a nonzero level-one modular form at
those weights by valence. Odd harmonic moments vanish by antipodality.
The fourth-moment identity gives horizontal perfection: a symmetric form
zero on all horizontal contacts has zero spherical mean square and is zero.

A symmetric form on the total 84-dimensional space vanishing on minimal
vectors has zero horizontal block by the preceding result. Evaluate it at
both (x,v) and (-x,v) for a mixed root contact. Averaging proves the vertical
form vanishes on all 24 D4 roots, so it is zero. Subtracting proves the cross
block applied to v annihilates all horizontal contact vectors for that color.
For the four simple roots, their contact vectors each span R80, so the cross
block vanishes. The entire lattice is thus perfect.

For positive eutaxy, use 80 literal horizontal contacts for **every** root
colour (remaining colours are in r*.out) with both signs, giving each
oriented vertical root the same total weight 1/12. Their outer products
sum to a block diagonal matrix with vertical part I4 and horizontal part S,
where S=(1/6) sum of the twelve per-colour mean horizontal outer products.
The goal is to represent I80−S using horizontal norm-8 contacts. This
condition is substantive: the mixed weights cannot simply be scaled down,
because their vertical part is already fixed.

The horizontal exact certificate exploits the native ten E8 blocks. The
240 doubled E8 roots split into fifteen orthogonal 8-frames up to sign.
For each frame, any pair of its roots on two distinct blocks, with
independent signs, is an actual BN80 vector of physical squared norm 8
(the replay verifies the matching lattice residues). Choose frame nine
as common orthonormal coordinates on all ten blocks. Every off-block
entry of I−S is supplied by the corresponding two-block norm-8 vector
with the correct sign; the resulting diagonal use in each of the two
blocks is the absolute off-block entry. What remains is ten symmetric
8×8 blocks T_i. Write u for a unit root direction (one representative
per sign pair). The local E8 root fourth moment gives
    T_i = sum_u ((10 uᵀ T_i u − tr(T_i))/30) u uᵀ .
The replay checks the fourth moment and the identity directly. The
minimum numerator before division by 30 is **35/96**, strictly positive.
The total coefficient over each one of the fifteen frames is exactly
tr(T_i)/15 by orthogonality of that frame.

It remains to realise these local diagonal terms by **actual two-block**
BN80 contacts; a single E8 root on one block need not be a lattice point.
For each frame distribute its ten block masses by symmetric edge weights
between distinct blocks, requiring the edge masses at a block to sum to
the given mass. This is possible whenever the largest mass is at most
half the total: partition a circle of circumference equal to the total
mass into ten consecutive intervals of those lengths; intersect each
with the translates of each other interval by a half turn. The result
is symmetric and has the specified margins and zero diagonal. Here the
balancing comparison is strict (trace slack 66563/1920 before division
by 15), so put a sufficiently small strictly positive weight on all
edges first and apply the same construction to the leftover masses.
At each block-pair edge distribute its weight over root directions
according to their independent within-block marginal proportions.
Use equally weighted plus and minus two-block contacts to cancel cross
terms and yield the desired diagonal contributions. This gives *strictly
positive* weight on every sign and direction in every native frame pair.

Thus the identity is a positive combination of mixed and horizontal
contact forms whose span is all symmetric forms. To see the latter with
precisely this support, pair signs give all horizontal off-block forms in
the common frame; horizontal direction squares are in the span of
pair sums between any three distinct blocks, and the local E8 root
directions are perfect by their fourth moment. The mixed forms span all
remaining coordinates by the preceding perfection argument. Every nonzero
symmetric infinitesimal form nonnegative on all contacts therefore has
strictly positive trace. Compactness of unit norm forms gives a uniformly
positive margin. The first-order Gram determinant expansion now proves
strict local optimality among all lattices, modulo scaling and isometry.

## Strict local optimum also among full-fibre periodic perturbations

Fix any finite-index subgroup P of D with F(P)=0 (replace an arbitrary
supercell by its common refinement with ker F). Consider deformations of
this finite full-H-fibre motif with fixed H, allowing independent phase
and vertical displacement of every site and an arbitrary deformation of
the vertical period basis, including horizontal components in its periods.
Let affine-linear drifts of the full finite-index periods be A vertically
and B horizontally; vertical and horizontal motif displacements u_v,h_v
are periodic mod P. Work infinitesimally after removing global translations
and vertical rotation (take A symmetric by polar decomposition).
At each root edge from v to v+r the two antipodal horizontal witnesses give
the inequalities
    r·(A r + u_{v+r}-u_v) >= |x·(B r + h_{v+r}-h_v)|
for every known contact x of that edge colour. They follow by differentiating
squared contact distances at min 8, with horizontal and vertical metrics
the Euclidean ones. Their left sides are nonnegative. Summing over all
sites and oriented roots cancels the periodic differences and leaves a
positive multiple of tr(A), since the D4 roots are eutactic.

If tr(A)=0, all the left sides vanish. Summing over sites for each fixed
root gives r·A r=0, hence A=0 by D4 perfection. The horizontal right sides
for four simple roots then force B r + h_{v+r}-h_v=0. Summing forces B=0,
and since simple roots generate D, h is constant. To see u constant,
take any nontrivial Fourier character chi of D/P. From the left sides
one obtains (chi(r)-1) r·u_hat(chi)=0 for each root. The roots with chi(r)≠1
must span R4: otherwise place them in a hyperplane spanned by three
independent roots (extend their span by roots), and all roots outside it
lie in ker chi. The finite replay verifies the roots outside each of all 24
such root hyperplanes contain an integral basis of D, contradicting the
nontriviality of chi. Thus u_hat(chi)=0, and all u are constant.

After quotienting translations and rotations, only the zero deformation
can have nonpositive volume derivative on the feasible tangent cone.
By the same compactness argument with first-order Taylor expansions, for
each fixed P this gives a strict local optimum for the finite motif
perturbations, measured by point density. Horizontal minimum 8 remains
attained in each fibre, and the stated exact periodic theorem fixes H and min8. It does not claim that
perturbations breaking individual H fibres, or allowing arbitrary changes
within H simultaneously with arbitrary finite period, are covered. There
is also no common quantitative radius asserted as P varies.

## Replay

Python 3 with NumPy and SymPy:
```
OPENBLAS_NUM_THREADS=1 python -B verify.py
OPENBLAS_NUM_THREADS=1 python -B eutaxy.py
```
Input `c84.json` is the parity chart, and `gram_and_block_map.json` is
the pinned BN80 Gram/atlas. Files c*.out and depth*.log are literal signed
contact residuals of denominator 64 in native coordinates, with Euclidean
metric raw-squared/2. The denominator and membership are checked directly
by exact basis inversion. The numerical determinant only proposes a basis
in the 24 hyperplane checks; each accepted basis has exact integer determinant
1 or -1. No weighted theta computation beyond the modular valence argument
is needed. Input decoding discovery (optional): `contacts.cpp` and
`input.txt`, cutoff 49152, 1400000000 nodes; e.g.
`./contacts input.txt phase2.in 49152 1400000000`. Cutoff/rank proof uses
only the literal witnesses.
