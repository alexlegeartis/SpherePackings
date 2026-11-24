# MW92: another large orbit excluded at invariant index 47⁶

At least **4,215,171,877 of 11,013,546,097** translation-invariant
index-\(47^6\) sublattices of Shioda's MW92 parent do not improve
the parent's center density. This is an increase of
**1,318,136,808** exclusions over the
[Round 21 screen](../../scoped_density_round21_20261124/92-mw-cusp-index6).
The proof still gives no conclusion about unhit profiles or general
sphere packings.

The additional row in `additional.json` lies at depth four in
the same six-jet quotient. Its squared norm is exactly 26, and its
normalized length-two ratio is \(19+24D\). Its verified normalizer
orbit comprises 276 previously unwitnessed length-two prefixes.
Those prefixes supply \(276\cdot47^4\) primitive profiles before
overlap with earlier deeper-prefix rows. The exact shortest-prefix
count after removing overlaps is:

| Witness depth | Shortest prefixes | Primitive profiles |
| ---: | ---: | ---: |
| 4 | 828 | 4,040,375,868 |
| 3 | 552 | 57,310,296 |
| 2 | 50,928 | 112,499,952 |

In addition, the earlier complete fifth-jet certificate excludes all
4,985,761 mixed profiles. The strict density comparison is
\(26^{46}<16^{46}\,47^6\).

## Complementary projection and an affine fiber

For discovery let \(M=D^4L\) and take an integral involutory Gram
isometry \(R\). The image \(M(I-R)\) is a lattice in the
negative eigenspace: a candidate \(w\) contributes \(w(I-R)\),
twice its orthogonal negative component. Row Hermite normal form
with quotient lifts gives a rank-46 basis for this image, plus a
rank-46 kernel over which each affine coset can vary. If a target
has squared norm at most 26, one of its two orthogonal components
has squared norm at most 13; using both signs of \(R\) therefore
requires projected squared norms at most 52 before halving.
Bounded sphere decoding in the projection and affine fiber
proposed the literal. Projection and decoding searches were capped;
neither is relied on as a complete short-shell classification.
The accepting criterion is the exact row norm and jet.

## Rational pole form

A depth-four norm-26 section has total intersection multiplicity
five with the zero section, four at the cusp at infinity. Consequently
there is one simple finite pole point \(a\in\mathbb F_{47^2}\);
after clearing the \(x,y\) denominators the section has form
\[
 X=P(t)/(t-a)^2,\qquad Y=Q(t)/(t-a)^3,\qquad
 Q^2-P^3=(t-a)^6(t^{47}-t),
 \quad \deg P=26,\ \deg Q=39.
\]
The finite pole can be reduced to either \(0\) or \(i\) using
affine transformations of \(t\) over \(\mathbb F_{47}\):
translation removes its real part, and a real nonzero scale
normalizes its imaginary coefficient when nonzero. All nonzero
scales in \(\mathbb F_{47}\) have sixth roots in
\(\mathbb F_{47^2}\), so the transformations lift to the curve.
This describes both pole charts, without saying either chart
is complete under any existing finite enumeration.

For the added witness the coefficients in `pole.json` give
\(a=40+34i\), with \(i^2=-1\). The independent pole replay checks
the polynomial identity, nonvanishing of both numerators at \(a\),
and all 92 exact height pairings against the pinned norm-16
basis sections. It does not rely on the finite-field reconstruction
used in discovery.

Replay `OPENBLAS_NUM_THREADS=1 python -B verify.py`, or
`OPENBLAS_NUM_THREADS=1 python -B verify_pole.py`
to check the extra geometric data as well. Optional explicit
arguments are Round 21's sixth-jet script, then the three
arguments documented by that script.
