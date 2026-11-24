# Split-shell designs for the 88D and 96D golden-cancellation branch

Use the hypotheses of the [local bound](README.md), in rank 88 or
96. In rank 96, **every nonempty signed shell with (q,s) fixed has a
spherical 3-design in each eigenspace separately** (normalizing by the
constant projected radius, and allowing multiplicity). The balanced
norm-12 shell is itself a spherical 3-design in 96 dimensions.
In rank 88 the projection to each eigenspace is in fact a spherical
5-design for every nonempty signed shell. These are necessary conditions,
not existence or nonexistence claims.

For bi-harmonic homogeneous polynomials of bidegree (a,b) on the two
eigenspaces, the weighted Hilbert theta series has weight (N/4+a,
N/4+b), with Fourier indices (n+j/sqrt(5))/2 where (n,j)=(q/2,s/2).
To check the full level-one transformation, upper translations in
the two embeddings follow from q even and <x,Ux> even, since
<x,Ux>=q(Ux)-q(x). Poisson summation in the two orthogonal summands
gives the inversion transformation, using Z-unimodularity of the
whole lattice. The product of the two Gaussian phases is trivial: its factor is
i to the power minus the total half-dimension (44 or 48), which
is divisible by four. Any sign from Fourier transformation of
the two harmonic polynomials also cancels, since their total
degree is even. The lattice need not be Galois-invariant.
Translations and inversion generate the Hilbert modular group
SL2(O_F), since O_F is norm-Euclidean: rounding the two coefficients
in the basis 1,U makes the absolute algebraic norm at most 3/4.
The Fourier coefficients of these weighted series are the separate
harmonic sums, and those at all n<4 and at n=4,j≠0 vanish.

[Aoki–Takemori, sections 2.2 and 3.1](https://www.kurims.kyoto-u.ac.jp/~kyodo/kokyuroku/contents/pdf/2100-07.pdf)
give the complete parallel ring and the mixed weight modules for
differences two and four. Their signed exchange is (-1)^k times
swapping the variables for parallel weight k, which explains why
the weight-5 generator is ordinary antisymmetric. Denote by G2,
G6 and G5 generators of the polynomial (signed-symmetric) ring R.
In the replay, G2 is the normalized weight-2 Eisenstein series,
G6 the weight-6 Eisenstein series with its constant removed by a
multiple of G2 cubed. This normalization does not affect the stated
Rankin–Cohen generating spaces (substitution of G6 plus a multiple
of G2 cubed makes the same module). The square of G5 is obtained
from the weight-10 Eisenstein series by canceling the constant and
first Fourier terms with G2^5 and G2^2 G6; the result is proportional
to G5^2 by the weight-10 basis. Recover the square root coefficients
in lexicographic order with leading coefficient 1 at (1,-1).
The parallel weight-15 antisymmetric (signed convention) generator
is obtained as the modular Jacobian of G2,G6,G5^2 divided by G5.
The Jacobian determinant has rows consisting of weight-multiple
of the form, derivative in the first embedding, and derivative in
the second embedding. It has weight 20 and is ordinary
antisymmetric; the ring gives its proportionality to G5 G15.
The n,j derivatives used in the code differ from the embedding
derivatives by a constant nonsingular change, irrelevant to spans.
Its computed leading coefficient is nonzero. The normalized
Eisenstein coefficients use the standard ideal-divisor formula
with multiplier 4/zeta_F(1-k) ([Dursthoff, Chapter 3](https://publications.rwth-aachen.de/record/685470/files/685470.pdf)).
The test of rational prime divisors is exact: 5 ramifies, 2 and 3
modulo 5 are inert, and in the split case the smaller valuation
is the valuation of gcd(n,(j-n)/2) in the basis 1,U.

For the modular dimensions the full parallel module is R + G15 R.
The difference-two module has generator degrees 7,8,11 and one
relation in degree 13. The difference-four module is free with
generator degrees 4,7. Therefore their weight-space dimensions
are computed by counting R monomials and, for difference two,
subtracting the dimension in degree shifted by 13.

Run python -B verifymod.py (standard library only). In rank 96 the
initial zero coefficients have full rank 9 for weight (24,26), so
this theta series vanishes identically, not just in a finite range.
For weight (25,25) the (n,j)=(4,0) and (6,0) functionals are already
in the span of the vanishing conditions. Thus both within-factor
and mixed covariance conditions hold on each of the two balanced
shells, and the equal squared projected radii give the same covariance
constant in the full space. Odd moments vanish by simultaneous
antipodality. In rank 88 the initial conditions have full rank at
weights for types (0,2),(0,4),(1,1),(1,3),(1,5),(2,4), and likewise
after exchange of the embeddings. At type (3,3) the balanced
norm-8 and norm-12 coefficients vanish as separate rows. No assertion
is made for arbitrary other bidegrees. Empty shells do not assert
the existence of a projected design.
