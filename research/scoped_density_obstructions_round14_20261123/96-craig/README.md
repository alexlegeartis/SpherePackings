# Rank-96 Craig phase and deformation exclusions

These are upper bounds on specific **families**, not a new density record or
an upper bound on all 96-dimensional packings.

Index coordinates by F_97. Let A be the integer augmentation lattice
(sum of coordinates zero), let M_j(x) = Σ_i i^j x_i modulo 97, and let
L_m be the sublattice of A with M_1, ..., M_(m−1) zero. The reference
metric is ordinary Euclidean. Vandermonde on 0,...,m−1 proves these
moments onto: solve with zeroth moment zero over F_97 and correct
integer coordinate 0 by a multiple of 97. Thus L_m has squared
covolume 97^(2m−1). For k distinct translates of period L_13 and
minimum squared distance q, squared center density is
k²(q/4)^96 / 97^25. Squared reference lengths in A are even.

**Claims.** Neither chart beats the existing [MW92–D4 96D lattice
construction](https://github.com/alexlegeartis/SpherePackings/blob/2724d234be553e7b7ffb6b6e730b7abde13e852b/research/mw92_cusp_d4_20261122/README.md), or a fortiori the
[stronger periodic construction](https://github.com/alexlegeartis/SpherePackings/blob/2724d234be553e7b7ffb6b6e730b7abde13e852b/research/mw92_cfr_taylor13_20261123/README.md):

1. Every union of translates of L_13 with phases in L_10/L_13,
   including arbitrary nonlinear and partial motifs.
2. Every affine subspace of L_9/L_13 as motif. In fact this also excludes
   an arbitrary such motif containing a full affine line with nonconstant
   M_9. A motif with constant M_9 is covered by 1 after translation.

Both claims hold for every positive-definite Gram deformation with
largest eigenvalue at most 26/25 after determinant-one normalization
and eigenvalue ratio at most 4/3. We do not exclude nonlinear L_9/L_13
motifs with varying M_9 that contain no such full affine line.

## Short lifting by moment coordinates

Under coordinate permutation i→ai+t the moments transform as
M_j' = Σ_(b≤j) binomial(j,b) a^b t^(j−b) M_b. All these permutations
preserve the reference length and each L_m.

If M_10 is nonzero, dilate to M_10=1 or 5 according as it is square
or nonsquare, and translate uniquely to M_11=0. Remaining dilation
freedom is ±1 and hence fixes M_12. The 194 normalized classes have
explicit signed-incidence differences of disjoint supports: q=20 for
77 of 97 square classes and 91 of 97 nonsquare classes, and q≤24 for
all 194. A recipe for these witnesses pairs split fibers of g(x²)
for degrees 5 and 6; the degree-6 pairs may differ in their x²
coefficient as well as their constant. Verification checks literal
moments, independent of the discovery recipe. At M_10=0, one q=22
vector with M_11 nonzero covers all such triples by dilation (11 is
coprime to 96) and translation in M_12. At M_10=M_11=0, the listed
q=24 witnesses cover every nonzero M_12 by dilation.

For a line direction with M_9 nonzero, translate to M_10=0 and divide
moments by M_9. Dilation acts on (U,V)=(M_11/M_9,M_12/M_9) by
(U,V)→(a²U,a³V). There are only 102 orbits, each with a literal q=18
signed difference. Scaling moments is allowed here because only
one nonzero residue on each line is needed; no scaling of the
short integer lift is performed. Every such full affine line
therefore yields a pair at squared reference distance at most 18.
The whole L_9/L_13 motif has at most 97^4 phases and loses at this
distance.

## A two-level bound for arbitrary motifs

Suppose reference minimum is at least 22. For each occupied first
moment u=M_10, define n_u as the maximum count of M_12 values with
M_10=u, M_11=v over choices of v. For two distinct occupied u,u'
and fixed second moments, normalizing the first two differences
takes the third difference by an affine bijection. Thus at most 20
third differences survive for u−u' square, or 6 for a nonsquare.
[Cauchy–Davenport](https://doi.org/10.1112/jlms/s1-10.37.30)
gives a difference set of size at least min(97,n_u+n_u'−1).
Consequently n_u+n_u'≤21 for square edges and ≤7 for nonsquare edges.

At most six occupied u have n_u≥4, because they form a clique in
the 97-vertex Paley graph, whose clique number the verifier computes
by normalizing an edge to 0,1 and searching its common neighborhood.
If max n_u≤6 the total of the n_u is at most 3·97+3·6=309.
If some n_u≥7, its 48 nonsquare neighbors are empty, and every
occupied n_u≤20 provided two first moments are occupied. Then the
total is at most 3·49+17·6=249. Thus in the multi-group case
k≤309·97=29,973; if only one first moment is occupied, k≤97².

For reference minimum at least 24, any occupied first moment has
only one second moment by the q22 lifts. Hence k≤309, or k≤97
if there is only one occupied first moment. At minimum at least 26
every nonzero phase difference is forbidden, so k=1. With no
distance constraint k≤97³. A q28 internal witness bounds the
singleton. These cases cover every possible finite motif.

For determinant-one G the first three bounds use the factor
λ_max(G) on squared reference lengths. For the singleton use the
affine orbit of the q28 vector below, whose mean squared length is
28 tr(G)/96. To bound this trace under condition ratio ≤4/3, scale
the least eigenvalue to 1. If M is their arithmetic mean, concavity
of log and log(4/3)≥2/7 imply geometric mean ≥exp(6(M−1)/7).
The arithmetic/geometric mean ratio is thus ≤(7/6)exp(−1/7),
at most 2401/2374 by bounding the exponential by its cubic Taylor
polynomial from below. Exact comparisons use squared lattice center
density (9/2)^96/(4·47^18), below the pinned periodic comparator.

## Intrinsic q28 orbit and small-index sublattices

The two supports in witness28.txt have difference in L_13 of squared
length 28. The affine group is 2-transitive: the mean of the orbit's
outer products is 28/96 times identity on the real augmentation.
The outer products also span all 4,656 symmetric directions.
Over F_389 use a primitive 97th root ω and nonzero Fourier coordinates
z_j=Σ_i x_i ω^(ij). Translation multiplies z_j z_k by ω^((j+k)t).
Taking the translation Fourier transform separates 97 blocks,
one per character j+k=s, each of dimension 48. Dilation reduces all
nonzero blocks to s=1. The verifier checks rank 48 exactly for
both blocks s=0,1 over F_389 with dilation rows 1,...,48 (the vector
is even under index negation). Thus the rational orbit spans, and
its minimum squared length strictly decreases to first order in
every nonzero trace-free Gram direction. This concerns the q28
orbit, not an assertion that L_13 has true minimum 28.

Every Euclidean index≤3 sublattice of L_13 contains a vector of
squared length ≤28. The three disjoint supports in tri28.txt have
equal moments through degree 12 and therefore form a q28 clique
after subtracting the first incidence vector. Index 1 and 2 follow
by pigeonhole. For index 3 extend the quotient map from L_13 to A
uniquely since the index of L_13 is a power of 97. If every cyclic
translate of the triangle exhausts the quotient F_3, then the sum
of its images is zero. The subtraction of the first incidence
has no effect mod 3. The sum of the three raw incidences belongs
to the mod-3 augmentation and its polynomial gcd with X^97−1
is X−1. Its cyclic shifts span the augmentation; thus the
quotient map vanishes, contradiction. These sublattices are
below the comparator since their covolume is at least that of L_13.

## Replay

Run python -B verify.py. All density tests are exact; printed decimal
ratios are illustrative. The largest Euclidean ratio to the pinned
lattice is 0.53183 (q28); the largest certified ratio under the stated
deformation class is 0.98357 (q22). Internal gaps, different periods,
and other carrier lattices are not excluded.

## Classification under translations and square dilations

There is only a one-parameter stretch rather than a general metric in
the translational and square-dilation-invariant rank-96 problem. Every
full-rank lattice in the real augmentation space invariant under
i→ai+t for square a and arbitrary t is, up to similarity and a real
isometry, a Craig L_m equipped with a determinant-one Legendre metric,
for m in 1,...,48.

To see this, put K=Q(zeta_97), O_K=Z[zeta_97], P=(1−zeta_97).
The translation T acting on rational points of the given lattice
has minimal polynomial Phi_97. Its fixed space under the square
dilations can be used to identify the rational representation with
the literal rational augmentation space: averaging over the dilation
group preserves rationality, and this fixed space spans the real
two-dimensional fixed space. Choose a rational fixed vector with
nonzero projections to both frequency blocks (each vanishing condition
is a proper subspace). Its translates give a cyclic generator.
The real intertwiner taking the standard orbit of e_0−(1/97)1 to
this orbit commutes with translations and square dilations. In
Fourier coordinates it is scalar on each of the two blocks, since
there are only square and nonsquare frequency orbits. Both scalars
are real because −1 is square, and signs can be absorbed by orthogonal
reflections. This intertwiner transports the original Euclidean
metric to a positive scale on each block.

Under the map x→Σ_i x_i zeta_97^i, translation acts by zeta_97
and square dilations by Gal(K/F) for F=Q(sqrt(97)), the fixed
quadratic field. Thus the transported lattice is a fractional
O_K-ideal stable under Gal(K/F). Away from 97 its valuations
are constant on sets above a common prime of F, and unramified
extension allows each such contribution to be pulled back from an
ideal over F. At 97 the unique prime P has relative ramification
index 48, so the ideal is a fractional ideal extended from F
times P^m. The class number of F is one: Minkowski's bound is
sqrt(97)/2<5, and elements (69+7sqrt(97))/2 and 10+sqrt(97)
have norms 2 and 3, respectively. Their conjugates generate the
other primes above 2 and 3 (the primes are distinct by
unramifiedness), so every ideal of norm at most 4 is principal.
Multiplication by the generator from F just rescales the two
frequency blocks independently.

The ideal P maps back to A, as seen from the differences e_i−e_j;
P^m maps back to L_m by expanding zeta_97^j at zeta_97=1.
For the range m≤48 this is exactly the first m moments condition,
using integral binomial coefficients and the triangular Newton change
of basis modulo 97. The Paley circulant J=(chi(i−j)) has
J^2=97I−11^T, and on augmentation is multiplication by the quadratic
Gauss sum in K, whose square is 97. Its ideal is P^48.
Multiplying by J therefore shifts m by 48 without changing shape
up to scale: J/sqrt(97) is orthogonal and commutes with the metric.
This gives precisely the stated range. The classification uses
standard cyclotomic ring and discriminant facts.

## Complete Legendre sweep for L_13

The two eigenspaces of J are ±sqrt(97), 48 dimensions each.
Write G_u for the metric with eigenvalues u and 1/u respectively.
Multiplication by 5604+569J is a lattice automorphism since its
inverse on augmentation is −5604+569J by the Pell identity
5604²−97·569²=−1. The automorphism identifies u with
(5604+569sqrt(97))²u; a nonsquare dilation exchanges u and 1/u.
It suffices to check 1≤u≤5604+569sqrt(97)<11209.

For a literal vector x, write q=x·x and g=x·Jx. Its squared
length is [q(u+1/u)+g(u−1/u)/sqrt(97)]/2. In frequency13.json
ten integral vectors of L_13 with g<0 cover consecutive intervals
from 1 to 11209: each supplied witness has squared length
less than 28.37 on its entire interval. For exact replay,
replace sqrt(97) with a slightly larger rational number in the
negative term to get an upper bound of the form au+b/u,
positive a,b. Convexity proves the bound from the two rational
endpoints. The resulting density is strictly below the fixed
96D lattice comparator for **every** L_13 Legendre parameter.

## Whole-frequency closure across 43 of the 48 classes

The same rational interval method applies to every moment class in
frequency_all.json, reusing a shared pool of exact short vectors.
The verifier checks lattice membership, exact Gram values, complete
interval coverage, and the strict density comparison for each
listed class. This closes the whole Legendre continuum for 43 of
the 48 Craig classes from the square-dilation classification. The
five classes not certified by this package are m=33, 34, 35, 37,
and 41. No claim is made for their uncovered parameter intervals.
