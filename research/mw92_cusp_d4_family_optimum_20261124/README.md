# Optimal fixed-MW92 cusp graph lattices; nonlinear cusp graphs reduce to 4D packing

These are rigorous **scoped** upper bounds, not a new upper bound for
arbitrary 96-dimensional packings. The arbitrary-coloring theorem misses
the current unrestricted 96D lower bound from the other construction
listed below, so this entire family cannot improve it.

Let M be the published Shioda rank-92 parent, using the same type-II
cusp and additive specialization rho as the existing MW92 + D4 graph.
Thus covol(M)=47^7, min²(M)=16, rho:M→F47² is onto, and its kernel K
has min²(K)=18 and covol(K)=47^9. For this result use also the checked
full-color witness: **each of the 2208 nonzero cosets of K has a
literal norm-16 representative**. The proof requires these literal
representatives, not any claim that a recorded stock exhausts the
norm-16 shell; the exact full-color verifier is a companion packet.

Let C0 be the center density of the known MW92 cusp + D4 graph,
namely (9/2)^48/(2·47^9).

## Arbitrary rank-4 linear graph tails

For any rank-4 Euclidean lattice T and any additive phi:T→F47²
(including non-surjective and zero maps), take the orthogonal graph
L={(x,y) in M⊕T: rho(x)=phi(y)}. The graph index is always 47²,
by surjectivity of rho, not by an assumption on phi. An overall
translation, or a constant affine offset in the graph equation,
does not affect the result; a constant offset is a translate by
a point of M.

Write q=min²(L) and V=covol(T). Then 0<q≤18 since K⊕0⊂L.
The pure tail T0=ker(phi) has index at most 47² in T and min²(T0)≥q.
By the exact 4D lattice Hermite bound attained by D4, any 4D
lattice with squared minimum m has covolume >=m²/2. Consequently
V≥q²/(2·47²).

Also any nonzero y∈T with phi(y)=0 must have ||y||²≥q. If its color
is nonzero, its literal norm-16 horizontal representative yields
||y||²≥q−16. Thus V≥(q−16)_+²/2 using the same Hermite bound.
Combining the volume inequalities gives the relative upper bound

    R(q) = 4 (q/18)^48 / max(q²/47²,(q−16)_+²).

On (0,q*] with q*=376/23 the first denominator wins, so R strictly
increases. On [q*,18], the logarithmic derivative is 48/q−2/(q−16);
it changes sign from negative to positive once at 384/23. Thus the
maximum is at q* or 18. Exact rational arithmetic in crossover.py gives
R(q*)<1/3, while R(18)=1. Since the existing D4 graph attains C0,
it is optimal in the whole stated family, not just among unscaled D4
tails tested at squared min 18.

## Nonlinear full-coset colored tails

Allow any finite-period arrangement of vertical points, each equipped
with a color c in F47² and carrying the complete horizontal coset
rho^−1(c). The arrangement need not be a lattice or a linear graph.
A color can occur multiple times, and when q≤16 different colors may
even share the same vertical point. For these arrangements q≤18 still
holds because every individual horizontal fiber is a full K coset.

Write nu for total tail-site intensity, counting colors separately,
and let Delta4 be the unconstrained occupied sphere packing density
supremum in R4, v4 the unit 4-ball volume. Each same-color subpacking
has squared separation at least q (matching horizontal points), and
there are at most 47² colors. If q>16, the whole vertical point set
has squared separation at least q−16, by the full-color literal
witnesses. Scaling these lower-dimensional packings gives

    nu <= 16 Delta4 / (v4 max(q²/47²,(q−16)_+²)).

The 96D center density is (q/4)^48 nu/47^9, so this family is bounded
by C0 times R(q)·Delta4/(π²/16), in particular C0·Delta4/(π²/16).

In fact the latter supremum is *attained as a supremum* at q=18.
To see the lower inequality, start with any 4D packing of squared
separation 2. Assign each vertical point a Voronoi cell of (2/3)D4
containing it, breaking ties arbitrarily. The D4 covering radius is
at most 1: round four coordinates to integers; if parity is odd,
change the coordinate of largest rounding error m≤1/2, giving
squared error ≤4m²+1−2m≤1. The scaled cell diameter is therefore
at most 4/3<sqrt(2), so each cell has at most one vertical point.
Color cells by the quotient modulo six times the cell-center lattice:
it has 6^4=1296≤2209 colors. Same-color centers have distance
at least 4sqrt(2), so same-color vertical sites have distance at
least 4sqrt(2)−4/3>sqrt(18). Distinct colors give horizontal
squared separation at least 16, and their vertical separation has
square at least 2. Lifting all sites to the corresponding full
horizontal cosets is thus a packing of exact squared minimum 18.

The lift can be made periodic at arbitrarily small density loss:
choose increasingly large finite near-optimal tail patches in boxes
of side 8m with a fixed-width guard around the boundary. Repeat
them by translations 8m·Z^4, a sublattice of 4D4, which preserves
the Voronoi colors. The fixed-width boundary loss tends to zero,
and the boxes can be chosen with tail density approaching Delta4
(for instance, clip increasing boxes of a periodic packing close
to Delta4). The upper and lower inequalities match as suprema.

## Numerics and exact scope

Table 6.1 of Cohn, de Laat and Salmon, *Three-point bounds for sphere
packing* (https://arxiv.org/pdf/2206.15373#page=30), gives the rigorous
rational upper bound Delta4≤0.636107333. Hence the 96D full-coset
family is at most 1.032 times the D4 graph (exact check in
crossover.py), below the unrestricted degree-13 covering route
(occupied density at least 6.250833919802144e−22, versus about
5.93710931559e−22 for D4 graph). The exact 4D lattice Hermite bound,
not a conjectured 4D unrestricted optimum, is used for the linear
case. See section 6.2 of the same paper for the distinction.
Arbitrary horizontal shifts not equal to these 2209 cosets,
different or deformed horizontal periods, and tilted metrics mixing
horizontal and vertical directions are outside the theorem.

## Reproduction and source identity

The parent and specialization map are the pinned published
[MW92 cusp + D4 graph](../mw92_cusp_d4_20261122/README.md), at the
cusp 18+2i (encoding 112). The exact Gram in gram.json.gz is extracted
from the source parent0.npz with SHA-256
`3c2a8629853cbd3f4ac37ad128a04b8a4c21db67a41888fed766ae840e0013f0`.
The included character matrix is the pinned evaluation character
of that same graph; the analytic minimum of the cusp kernel and
its norm-18 certificate are inherited from the parent construction.
The audit uses one integer row per color and no external shell search:

    python -B verify_colours.py
    python -B crossover.py

The color table covers all 2208 nonzero values by literal q=16 rows.
In particular the proposed scale 1/3 tail with color matrix
[[9,39],[17,3],[22,9],[23,36]] has a matched vector over the tail
(-1,-1,0,0), color (21,5), of squared norm 16+2/3,
not 18. Its conditional ninefold gain is therefore false. The
published lattice and unrestricted density records remain unchanged.
