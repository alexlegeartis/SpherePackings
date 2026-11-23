# Rank 96: exact Craig exclusions beyond the native scalar charts

This is a scoped exclusion, not a new packing. On the rank-96 augmentation
lattice A of Z^97 use L_m = {x in A: sum i^j x_i = 0 mod 97 for j<m}.
The squared covolume is 97^(2m-1). Against the live MW92 glued 96D
periodic packing use squared center density 4^96/(47^14 D), with
D = 10950866250772257481/50625000000000000. All comparisons in verify.py
are rational.

## Small indices and nonlinear phases

For h|96 and m<=h, take the r=96/h nonzero fibers F_j of X^h in F97.
The h-point fiber powers agree through degree h−1. Differences of
incidences are in L_m and have squared norm 2h. Subtracting one fixed
incidence gives a clique P of r points in L_m. Any index k<r sublattice
contains a difference from P. For k=r, assume no such short vector at
all. Let G=L_m/L' and extend the quotient map uniquely to A: the index
[A:L_m] is a power of 97, invertible on G because |G|=r divides 96.
On each cyclic coordinate translate of P the images must exhaust G,
so their sum is constant. Their raw sum at one translate is
1−e_t−r*1_(t+F_0), an augmentation vector. Taking consecutive
differences and using rG=0 forces all e_t−e_(t+1) to have zero image,
contradiction. Thus k<=r (inclusive) gives squared min<=2h.
In particular no index <=8 Craig L11 or L12 sublattice beats the live
bound in the Euclidean metric. For L11 index 1 use a q22 witness, index
2 use the odd-length cyclic orbit of this vector (sum zero), and
indices 3..8 use q24. L8 indices up to 12 also lose. A literal q30
vector in L15 lies also in L14 and shows both raw parents lose.

For nonlinear phase unions in L11 with period L13, **all 9408 nonzero
classes** of L11/L13 have q22 or q24 literal witnesses. For the 9312
classes with M11 nonzero, one q22 representative has nonzero M11 and
dilation covers all nonzero M11 because gcd(11,96)=1. Translation
changes M12 by 12*t*M11 and realizes all M12. For the remaining 96
classes use q24 vectors in L12: the replay tests their M12 values after
dilation and finds all 96 values. Any construction with this fixed
period and minimum at least q26 has at most one phase; it cannot
improve at q26 or q28. This doesn't exclude a better period, genuine
shifts outside L11, or q30 or higher for L13 itself.

## Frequency metric

The Legendre circulant K=(chi(i-j)) obeys K²=97I−11^T, hence on A over
R it has eigenvalues ±sqrt(97), each 48 times by trace. The positive
determinant-one metric G_u = ((u+1/u)I+(u−1/u)K/sqrt(97))/2 scales the
two eigenspaces by u and 1/u. K takes L_m into L_m by convolution and
the binomial theorem. Multiplication by 5604+569K is an automorphism:
inverse −5604+569K on the augmentation lattice by Pell norm −1. The
metrics at u and (5604+569sqrt(97))²u are therefore isometric. A
nonsquare dilation exchanges u and 1/u. The full analysis reduces
exactly to 1<=u<=5604+569sqrt(97).

Literal contact witnesses certify no improvement on rational compact
subintervals starting at u=1: upper endpoints are 2.49,2.71,1.87,1.84,4
for m=8,9,10,11,12. For each, choose the sign of x^T K x nonpositive via
nonsquare dilation and use sqrt(97)<197/20. The resulting rational
upper estimate is convex as a function of u, so endpoint checks suffice.
The exploratory floating point LLL scan did not certify or exclude any
remaining parameter. It is included as a search log, not as evidence of
a lower bound or a general arithmetic obstruction.

The auxiliary exhaustive equivariance scan compares full fibers of
X^(m mod d)g(X^d), monic g. No m11 pair found in d=3,4,6,8; no m13
pair found in d=4,6,8,12, nor among 1533939 normalized odd polynomials
with fully split zero fiber tested separately. The m15,d6 scan did
find a valid pair, checked independently by verify.py. The q22
triangle screen over sampled stocks is not an exhaustive shell
certificate. No compute remains active.

Archived replay SHA-256: `721e6e159fc8458c55a06c610a808f4052031930dfb4309acf5777aa1893d1f7`.
