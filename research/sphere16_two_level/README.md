# Large exact metric chambers around extremal Type II lattices

Let L be an extremal even unimodular lattice of rank n=48 or 72, of
squared minimum m=6 or 8. Its original metric maximizes density,
strictly up to common rescaling, in the following **scoped families**:

* Every positive definite metric with at most **three distinct eigenvalues**,
  with arbitrary multiplicities and eigendirections, whenever their ratio is at most
  **15/4 in 48D** or **3 in 72D**. For two levels the rank-dependent
  thresholds in the table are stronger: a rank-one stretch allows ratio at least **6**
  in both ranks, and a rank-one contraction allows **30,000** and
  **38,000**, respectively.
* Arbitrary contractions in up to k directions. There is **no
  requirement that these k contracted eigenvalues agree**. Sample
  certified condition ratios, rounded down from the exact table:

  | k | 48D | 72D |
  |:--:|:---:|:---:|
  | 2 | 371 | 415 |
  | 3 | 82 | 88 |
  | 4 | 36 | 39 |
  | 5 | 20 | 21 |
  | 6 | 12.7 | 13.6 |
  | 8 | 6.22 | 6.92 |
  | 10 | 3.60 | 4.10 |
  | 12 | 2.24 | 2.68 |

* Arbitrary positive semidefinite stretches in up to k directions
  with possibly unequal eigenvalues. Exact ratios:

  | k | 48D | 72D |
  |:--:|:---:|:---:|
  | 2 | 5.5 | 5.2 |
  | 3 | 5 | 4.75 |
  | 4 | 4.7 | 4.4 |
  | 5 | 4.5 | 4.1 |
  | 6 | 4.3 | 3.9 |
  | 7 | 4.1 | 3.75 |

This enlarges the exact shape controls at high condition ratios, not the
condition-number balls for unrestricted spectra. It neither constructs
a new lattice packing nor bounds arbitrary packings.

## Shortest-shell projection tests

In ranks 48 and 72 the extremal shortest shell is a spherical
11-design. To see this, the degree-j harmonic theta series has weight
n/2+j and vanishes to order at least m/2; division by Delta^(m/2)
therefore leaves a holomorphic level-one modular form of negative
weight j-12 for even 2 <= j <= 10. Odd moments vanish by antipodality.

For any orthogonal projection P of rank r, let X=||Pv||^2 for a
unit direction v in the original shortest shell. Its first five
moments m_j equal those of Beta(r/2,(n-r)/2): this follows from the
squared lengths of the two orthogonal projections of a standard
Gaussian and the independence of direction and total length. Thus
m_0=1 and m_(j+1)/m_j=(r+2j)/(n+2j).

Put H=(m_(i+j)), J=(m_(i+j+1)) for i,j=0,1,2. The rational witness
u in the table is chosen with det(J-uH)<0. In particular J-uH is
not positive semidefinite. Since its quadratic form on polynomial
coefficients is the shell average of (X-u) times a quadratic square,
some shortest direction has X<u. This applies for every P, without
exhausting any eigendirections. Applying the same argument to I-P
also gives a shortest direction with X>1-u_(n-r).

## Two eigenvalues

Normalize the metric to G=I+(Q-1)P, where Q>1 is the condition ratio.
Its squared minimum is less than m(1+(Q-1)u_r) by the projection
test. The determinant root is Q^(r/n). The exact comparison at
each table endpoint Q_max is

    (1+(Q_max-1)u_r)^n < Q_max^r.

By strict concavity of Q^(r/n), comparison with the same affine
function holds for the whole interval 1<Q<=Q_max. Density cannot
increase: relative to the original, it is the n/2 power of the quotient
of the squared-minimum ratio by the n-th root of the metric determinant.

## Contractions of bounded rank, unequal eigenvalues allowed

Normalize G=I-A where A is positive semidefinite of rank at most k<n,
with nonincreasing eigenvalues alpha_1,...,alpha_k (pad by zeros).
Set T=1-1/Q_max and suppose alpha_1<=T. For j=1,...,k apply the
projection test to the span of the top j eigendirections of A, using
c_j=1-u_(n-j). As A is at least alpha_j times that projector in the
positive semidefinite ordering, the squared minimum is at most
m(1-max_j alpha_j c_j). Put h=(max_j alpha_j c_j)/T. Then

    alpha_j <= T min(1,h/c_j),
    (det G)^(1/n) >= D(h) = product_j(1-T min(1,h/c_j))^(1/n).

Here 0<=h<=c_k, and c_1<...<c_k. On each of the intervals cut out by
0,c_1,...,c_k, the function D is the geometric mean of n positive
affine functions (the additional n-k functions equal 1). It is
concave on each interval, by concavity of geometric mean. Thus the
integer/rational inequalities D(c_j)>1-Tc_j at all the positive
endpoints, and equality at h=0, give D(h)>1-Th for every h>0.
This is precisely the needed strict determinant comparison for
general contractions, not an interpolation of density between
different lattices.

## Three distinct eigenvalues and arbitrary low-rank stretches

For three distinct eigenvalues normalize them to 1, 1+(Q-1)t, Q with
0<t<1. Let the high and middle multiplicities be r,s. The projected
moment variable is Y = U+tV where U and V are the squared projections
of the unit shell direction to the orthogonal high and middle spaces.
Their mixed moments are

    E U^i V^j = (r/2)_i (s/2)_j / (n/2)_(i+j),

by the same Gaussian direction calculation. The Hankel test for Y
therefore applies for the same reason as the Beta test. On a t-interval
[a,b], any pointwise Y at t is bounded above by Y at b, and the metric
determinant is at least Q^r(1+(Q-1)a)^s. Select the largest rational u
of denominator 10^9 with (1+(Q-1)u)^n smaller than this determinant;
use the moments at b and certify that the matrix of (Y-u) times a
quadratic square is not PSD, by a negative principal minor. If this
fails, subdivide the interval until all boxes close. Because the
projection comparison has an actual shortest direction as witness,
and the geometric mean of the eigenvalues is concave in Q on the
fixed shape, proving the comparison at the desired Q_max proves it for
all 1<Q<=Q_max, including those nearer the scalar endpoint.

For a k-direction stretch, normalize the metric to I+(Q-1)A, with
the ordered nonzero eigenvalues of A in [0,1], largest 1. Pad by
zeros if the rank is smaller than k. For each box of the other
k-1 eigenvalues, compute the same moment matrix at their simultaneous
upper endpoints and use the metric determinant of their lower
endpoints. If the box intersects the ordered domain, tightening its
interval endpoints by the order constraints is valid. The moment
variable is the linear combination of the k squared projections onto
the corresponding orthonormal axes; its mixed moments are again the
Gaussian projection formula, now with each factor of rank 1.
The verifier subdivides boxes, tests exact signs (no solver or
floating point), and uses the same concavity for intermediate scales.

## Replay and reference

Run python -B verify.py with Python 3. No floating point or
third-party dependencies enter the verification. In inputs.py,
the two integers in each pair are the numerators of u_r (denominator
10^12) and certified Q_max (denominator 1000), in ascending r.
The contraction list has one Q_max for each k=1,...,12 (denominator
1000). The exact test checks all determinants and both families of
density inequalities. Run python -B three.py for the three-level
claims (55,799 and 134,153 exact interval nodes), and
python -B stretch.py for the unequal low-rank stretch claims; an
optional argument 48 or 72 selects one rank.

The design argument uses the standard weighted-theta theorem; see
[G. Nebe, Boris Venkov's Theory of Lattices and Spherical Designs](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf).
See the [earlier unrestricted-spectrum shape chambers](../p48_gamma72_metric_and_equality_20261123/README.md)
for the different scope Q<=2 and Q<=19/10.
