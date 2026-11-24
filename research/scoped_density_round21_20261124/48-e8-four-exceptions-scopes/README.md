# Four exceptional E8 slots: exact determinant subregions

Consider the positive-definite 48-dimensional E8 seven-slot zero-sum
half-period model with slot weights (a,b,c,d,1,1,1), where a<0 and
b,c,d>=2. For arbitrary half-period motifs, center density is at
most 16,800 if any of the following holds:

* D<=1/3500;
* D>=128/1875;
* all three heavy weights lie in [2,5/2].

Here D=bcd+a(bc+bd+cd+3bcd) is the scalar determinant. The proof
does not cover the outstanding interval 1/3500<D<128/1875
outside the bounded heavy-weight cube.

## Normalization and proof

The covolume of the doubled period in the scaled metric is 2^24 D^4.
Here g is squared minimum, N<=2^48 counts occupied half-period
phases, and the radius is sqrt(g)/2.

Along one common E8 root, restrict to integer vectors of slot
coefficients (-(u+v+w+3t),u,v,w,t,t,t). The scalar form before
multiplication by 4 has determinant 3D (matrix
diag(b,c,d,3)+a(1,1,1,3)(1,1,1,3)^T). The exact rank-four Hermite
constant sqrt(2) implies g^4<=3072D. Before accounting for the half-period direction, this gives density
N(3/4)^6 D^2. The lower cutoff uses the improvement below. In each determinant
slab it also gives the factor bound (3/4)^4 (g/8)^8 for density
divided by N. A single opposite root pair on a and weight 1 gives
g<=4(1+a)<4, proving the uppermost cutoff D>=6.

At fixed D, a=(D-bcd)/(bc+bd+cd+3bcd) separately decreases in
b,c,d. For example its b partial derivative has numerator
-[(cd)^2 + D(c+d+3cd)] and denominator (bc+bd+cd+3bcd)^2.
It is separately convex as well. Thus estimates involving only a
and weight-1 slots can be evaluated conservatively at b=c=d=2.
These comprise short scalar periods, the two-slot argument for
a+1, and the exact rank-three Fourier quotient with changing slots
a,1,1,1. The method is the same as in
[the completed three-exception packet](../48-e8-three-exceptions-complete) with shape41 now unavailable.
All Fourier coefficients and rational bounds are checked in replay.

Finally on the bounded cube use the root direction
(-9,1,1,1,2,2,2). Its squared length 4(81a+b+c+d+12) has a maximum
at a cube corner for fixed D by separate convexity, and increases
with D. Rational comparisons on twelve slabs close the seam.

Replay: python3 -B replay.py (requires a standard C++20 compiler).

## Primitive-direction improvement

Any primitive integer slot coefficient vector (with total zero)
defines a copy of E8 in the half-period carrier. If its doubled-period
root squared length is p, nonzero parity cosets on this E8 direction
have representative squared lengths at most p/4 for anisotropic and
p/2 for isotropic labels. Primitivity ensures that the E8/2E8 label
subgroup embeds in the phase group. Thus, when g>p/4, the number of
centers in each directional fibre is at most 16, by the maximal
dimension four of a totally singular subspace; when g>p/2, each
fibre is a singleton. Also g<=p. For any exponent r>=8 the three
ranges therefore imply N g^r <= 2^40 p^r. This applies to the
rank-four Hermite shortest vector (primitive) for r=24, giving
density <=2^40 (3/4)^6 D^2. It also applies to the fixed cube period,
which is visibly primitive, for each exponent 8 or 24 arising
in the elementary density envelopes. This yields the improved
D and cube cutoffs stated above.
