# Periods in the all-good E6×E6 packet mask

We study the rank-five quadratic quotient W=E6/(3E6+R) appearing in
the 12D equal-scale ternary mask (equivalently (E6/3E6)/R, where R
means the mod-three radical). An all-good code is a subset of W×W
whose two projection difference types for every distinct pair are
(1,3), (3,1), (2,2), (2,3), (3,2), or (3,3). Type 0 is the zero
vector, 1 and 2 are quadratic values, and 3 is nonzero isotropic.
Each code label can carry a full nine-site radical packet.

**Results.** Every all-good code with a nonzero translation period
has at most 93 labels. If its translation stabilizer has dimension
two then it has at most 63 labels, and the literal supplied witness
attains 63 (567 motif sites). The unrestricted 103-label target
remains open and any such code must have no nonzero translation period.

## First reduction: the period is isotropic in both projections

Suppose the code has a nonzero period (u,v); its cyclic group is a
line L. Its type must itself be good. Every allowed difference of
L-cosets must admit all three translates as good differences
(unless the cosets coincide). The exact Delsarte dual in
`period1.json` gives, for period types (1,3), (2,2), (2,3),
respectively a cap of 9, 9, 15 cosets, i.e. 27, 27, 45 labels.
Swapping factors covers the reversed types.

The spatial classes in the computation use the cyclic tuple of the
three types of x+i*u and y+i*v for i=0,1,2. In each factor those
types determine the Gram and dependency data for the ordered pair
(x,u): the zero type recognizes collinearity, and the quadratic
values determine the polar pairing. Witt extension therefore
identifies each fiber of this tuple with a stabilizer orbit. The
common cyclic quotient shift accounts for cyclic rotation. Ordinary
real character inequalities on the quotient are weighted by the
provided rational nonnegative dual. The verifier needs only a
subset of all possible quotient character inequalities; it computes
them as 3 times the number of zero residues minus the number of
quotient group elements, without any Fourier rounding.

Hence any stabilizer of a code with more than 45 labels must have
all its nonzero vectors isotropic in both projections. Its
projections are injective, and the Witt index of W is two, so
the stabilizer has dimension at most two.

## Isotropic lines and planes: covector partition

Align the two projections of an isotropic translation subspace U
by Witt isometries, so the period is diagonal. In each W-factor
use Witt coordinates (a,b,X), where U shifts a, b is its dual
covector, and Q(a,b,X)=a·b+Q_res(X). For dim U=1, Q_res on
a rank-three space has the model xy-z²; for dim U=2 the residual
is -z². The anisotropic square class is determined by the rank-five
shell counts (1,72,90,80), after removing two hyperbolic planes.

Between cosets, the pair of quadratic values as the period shift
varies is an affine image of U. In addition to never vanishing
as an actual vector in either factor, the pair of quadratic values
must lie in {(0,0),(0,1),(0,2),(1,0),(2,0),(2,2)}. The only
nonconstant affine lines in this six-value table are the two
axes and the antidiagonal a+b=1. Affine rank two is impossible.
Thus for the pair of covector differences, either one covector
is zero, or the two are nonzero and opposite. In the opposite
case Q_total of the difference is 1. In an axis case, the
zero-covector factor has quadratic value zero.

## Plane periods: sharp cap 63

For dim U=2, an axis case is impossible: a zero dual covector and
zero quadratic value force zero residual in that factor, so
a period shift makes its actual vector zero. Thus the sum of the
two covectors is fixed, and after translating it can be taken
as zero. Across unequal covectors the total quadratic distance
is 1. At equal covectors, the two residual differences must
both be nonzero (otherwise a period shift kills one projection),
so both have value 2 and the sum is again 1.

Write the two W-factors as (a,b,t),(c,-b,s). The quotient of
b+(-b)=0 by the diagonal period has the six coordinates
(a-c,b,t,s), and the induced quadratic form is
(a-c)·b-t²-s². All distinct occupied cosets have quadratic
distance 1. Translated to include zero, the polar Gram of n
remaining vectors has diagonal 2 and off-diagonal 1, hence is
I+J over F3. Its rank is at least n-1. It cannot have n>=8,
and for n=7 its rank is 7 because n+1!=0 mod 3. Thus n<=6:
there are at most seven cosets, i.e. 63 labels. The supplied
witness has exact stabilizer size 9 and 63 labels.

## Line periods: cap 93

For dim U=1 the covector-pair fibers lie in F3². Two distinct
fibers must differ by an axis direction or the antidiagonal
direction. This graph has clique number three. If a clique of
fibers is not all on one antidiagonal and contains at least
two fibers, every fiber is incident to an axis edge. These
facts are checked by direct finite enumeration.

Within each fiber the residual projection to each three-dimensional
factor is injective: if one residual difference vanishes, a
period shift kills the corresponding entire projection. The
residual pairs form a rank-three product good code, with exact
two-point cap 18. In an axis pair, the zero-covector residual
factor has all cross differences nonzero isotropic. Fixing a
single coset in the opposite fiber bounds the first fiber by
the eight nonzero isotropic vectors in that projected factor.
This applies to both endpoints. Thus axis-incident cases
have at most 24 cosets; a single fiber has at most 18.

It remains that two or three fibers lie on one antidiagonal.
In coordinates (a,b,X),(c,-b,Y), across distinct values of b
the total quadratic distance is constantly 1. Let S_b be
the linear span of all differences of residual pairs (X,Y)
in fiber b, inside the nondegenerate residual direct sum of
dimension six. Double subtraction of the constant cross
quadratic distances gives S_b perpendicular to S_c for b!=c:
the discarded term (a-c)·b has zero double difference
because the covector is fixed within each fiber.

Writing r_b=dim S_b, a fiber has at most min(18,3^r_b)
cosets. At r_b=3 this improves to 13: more than half of a
27-point affine space would realize every difference in
its span, forcing a full rank-three additive good code in
the rank-three product with both projections bijective.
The latter is impossible by the exact small quadratic-form
replay. If S_b is totally isotropic, the fiber has at most 3:
its allowed total values are then (0,0) only, so its first
residual projection spans a totally isotropic affine line.

If some r_b>=5, every other rank is <=1, giving at most 24
cosets. If the maximum is 4, every other rank is <=2.
If two others have rank 2, they both equal the same
two-dimensional orthogonal complement of the rank-four span,
and mutual orthogonality makes that complement totally
isotropic, giving each such fiber at most 3. Otherwise
the bound is 18+9+3=30. If exactly one span has rank 3
and all others rank <=2, the bound is 13+9+9=31.
If two spans have rank 3, each is the full orthogonal
complement of the other, and the remaining span is in
their totally isotropic intersection: the bound is
13+13+3=29. Finally when all ranks <=2 the bound is 27.
Thus there are at most 31 cosets, i.e. 93 labels.

Replay these results with Python and NumPy:
```sh
python -B verify_period1.py
python -B verify_small.py
```

## Occupancy five absent, A=110–117

Retain the previous equal-mask proof that any record-size
submotif using only packets of size at least four has
103<=A<=117 occupied packets. If no packet has size five,
let J count packets of size at least six. The number of sites
is at most 4A+5J, so a record needs J>=ceil((923-4A)/5).
Each of those J packets is universal-good in this submotif.
The new independently exact marked certificates exclude
A=110,111,112,113,114,115,116,117 at their respective J
thresholds 97,96,95,95,94,93,92,91. This leaves A<=109
in the zero-five branch. Use the exact preceding round
marked model and verifier with the additional certificates:
```sh
for p in mcert*.json; do python -B verify_marked.py "$p" || exit 1; done
```

### Remaining zero-five cases reduce nearly entirely to the all-good target

The singleton-bad graph on the size-four packets is triangle-free,
by the canonical forbidden-offset balance on a triangle: after
translation three four-sets would have to be pairwise disjoint
in a group of size nine. All J higher-occupancy labels are good
to these packets and to one another. A triangle-free graph with
at least 6 vertices has an independent triple (the elementary
R(3,3) bound); with at least 9 it has an independent four-set.
For the second bound, a counterexample on nine vertices has every
degree at most three, as every neighborhood is independent; every
degree is also at least three, since otherwise at least six
nonneighbors yield an independent triple, which can be joined to
the vertex. The total of nine degrees equal to three is odd.

Check each J>=ceil((923−4A)/5) for 103<=A<=109. Adding the guaranteed
independent set of fours to the J universal-good labels always
gives a 103-label all-good subcode, except for:
- A=107, J=99 with independent number three on its eight fours.
  All 99 high packets must be full, since 4*107+5*99=923.
- A=109, J=98 with independent number four on its eleven fours.
  At least 95 of the 98 high packets must be full, since
  4*109+5*98=926 allows at most three missing sites.

A further phase consequence in either exceptional case is that
the triangle-free graph has a short odd cycle with nonzero total
gain for the unique forbidden offsets. In the eight-vertex case
there must be a five-cycle: a shortest odd cycle of length seven
would have no chords, and the remaining vertex has at most two
neighbors on it (otherwise two neighbor cyclic distances create a
triangle or a five-cycle). Deleting those neighbors from the cycle
leaves an independent triple and the remaining vertex can be joined
to it. In the eleven-vertex case there must be a five- or seven-cycle:
otherwise a shortest odd cycle is a chordless nine- or eleven-cycle.
An eleven-cycle already has five independent vertices. A vertex
outside a nine-cycle has at most two neighbors on it (the only
allowed cyclic distance between two neighbors is two), and avoiding
those neighbors again gives an independent four-set on the cycle
which can be joined to that vertex.

Finally a zero-gain odd cycle of length five or seven is impossible
for fours. Zero gain permits aligning all packet phases so that
each edge requires disjointness of its two four-sets. Along any
two-edge path the endpoints then have at least three points in common
(their union avoids the middle four-set in a plane of size nine).
After four or six edges at most two or three points of the initial
four-set have been lost, so the closing edge cannot be disjoint.
This is a necessary cycle condition, not a dismissal of either
exceptional branch.
