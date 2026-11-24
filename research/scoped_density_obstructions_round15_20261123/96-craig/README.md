# Rank 96: a robust exclusion around all square-affine Craig shapes

In the 96-dimensional augmentation space indexed by F_97, put
`L_m={v in Z^97 : sum v_i=0, sum i^k v_i=0 mod 97 for 1<=k<m}`.
Let J be the quadratic-character circulant; G_u is the determinant-one metric
with eigenvalues u and 1/u on its two 48-dimensional eigenspaces.

**Results.** For all m=1,...,48 and all u>0, the center density of (L_m,G_u)
is less than 543/1000 = 54.3% of the
[MW92–D4 rank-96 lattice center density](https://github.com/alexlegeartis/SpherePackings/blob/2724d234be553e7b7ffb6b6e730b7abde13e852b/research/mw92_cusp_d4_20261122/README.md).
The [prior classification](https://github.com/alexlegeartis/SpherePackings/blob/d44afaa088b1c3b64f22dbb88e767b1d10cb925b/research/scoped_density_obstructions_round14_20261123/96-craig/README.md)
shows this exhausts lattices invariant under the standard permutations
`i -> ai+t` with a square, up to similarities and isometries.
Moreover *arbitrary*, potentially symmetry-breaking, metric deformations of
these lattices cannot beat MW92–D4 when the generalized eigenvalue ratio
relative to any G_u is at most 11/8. No claim is made about other lattices
or unrestricted packings. In the deformed neighborhood the certified bound is
99.8% of the MW92–D4 comparator.

## Finite exact envelope proof

The supplied vectors have integral coordinates, and membership is checked
through moment m−1. For any nonzero witness v write q=v·v and g=v·Jv≤0.
Its squared norm in G_u, u≥1, is at most au+b/u, where
a=(q+g/S)/2, b=(q−g/S)/2 and S=4,924,428,901/500,000,000 > sqrt(97).
The affine permutations preserve each L_m, and on the augmentation space
the Pell unit 5604+569J has inverse −5604+569J, both preserving L_m
as integral cyclotomic multipliers. This identifies u with
`(5604+569 sqrt(97))² u`. A nonsquare dilation exchanges u and 1/u.
Consequently [1,11209] covers the full stretch domain.

Put w=u². The minimum of the witness bounds is the lower envelope of
the affine functions aw+b, divided by sqrt(w). The breakpoint w between
two consecutive lines is rational. On each segment au+b/u is convex in u,
so a bound on the fourth power of the distance is obtained from the two
exact endpoint tests (aw+b)²/w. A bound U for that fourth power gives
a bound (U/16)^48/97^(2m−1) on squared center density; this is compared
to (543/1000)² times the squared MW92–D4 comparator. The strict comparison
holds on all intervals for each of the 48 classes.

## Arbitrary deformation by averaging

Let H be a positive-definite metric on L_m. Average H by conjugating with
the square-affine coordinate permutations, obtaining `bar(H)`.
Because those permutations preserve L_m, the original minimum is bounded by
the averaged metric minimum: take a shortest vector for bar(H), and one
member of its orbit has H-length at most its orbit-average squared length.
The averaged metric is a scalar times some G_v by the two Fourier orbits;
the scalar cancels when computing its normalized density. Thus, with
D the MW92–D4 comparator,
`density(L_m,H) < (543/1000)*sqrt(det(bar(H))/det(H))*D`.

Suppose now that the eigenvalues of `M = G_u^(-1/2) H G_u^(-1/2)`
lie in [1,11/8] after a common scaling. Averaging commutes with conjugation
by the G_u square root. By the arithmetic-geometric mean inequality,
det(bar(H))/det(H) is at most (arithmetic mean / geometric mean)^96
of the eigenvalues of M. To bound that ratio, the chord below log on
[1,11/8] gives mean-log ≥ c(MEAN−1), with c=(8/3) log(11/8).
Maximizing log(MEAN)−c(MEAN−1) yields c−1−log c.
Writing x=3/8, the alternating sum c0=Σ_(j=0,...,9) (−1)^j x^j/(j+1)
satisfies c≥c0 and c≤1. With t=1−c0, the power series bounds
c−1−log c ≤ t²/2+t³/3+t⁴/4+t⁵/[5(1−t)].
Forty-eight times this last expression is less than the alternating
20-term lower bound for log(998/543), using y=455/543 as argument of
log(1+y). The replay verifies the rational comparison, so the determinant
factor is less than 998/543.
The deformation need not commute with any coordinate permutations.

## Replay and discovery

Run `python -B verify.py` using only the standard library. This checks all
moments, scalar products, rational coverage, density comparisons and the
rational logarithmic inequality. The vectors were located in rational
isotypic subspaces, using even/odd subspaces with progressive metric
reduction near the Pell endpoint, and signed shell search in L_7.
No numerical lattice reduction enters the proof.
