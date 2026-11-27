# Addendum: exact 26D ceiling over all X13 dual phases

Take the period \(X_{13}\) and phase alphabet \(X_{13}^*/X_{13}\)
defined in the companion README. This includes all \(3^{15}=14{,}348{,}907\)
phase labels as possible centers. For an arbitrary subset of them
and any radius, the center density is at most \(1/\sqrt3\).
The 2187-phase witness in the README attains it. This excludes
new records in this fixed isotropic alphabet, including all nonlinear
ternary motifs, not arbitrary real shifts or a different period.

## Unit radius (\(\mu=4\))

The exact metric makes each occupied bank a ternary code of
length 13, minimum Hamming distance 6, and with one color per
coarse word. Every bank has at most 1449 words by the published
[Gijswijt–Schrijver–Tanaka \(A_3(13,6)\) bound, Table 1](https://ir.cwi.nl/pub/13374/13374D.pdf).

Across any two banks only bracket 2 is allowed, because the
three consistent squared torus distances are \(26/9,32/9,38/9\).
The affine equation of the README at each cross edge may have
a different constant because of the wrap. This does not change
the rectangle difference or the cancellation of each color
between incoming and outgoing edges. Thus the difference spans
\(U_0,U_1,U_2\subseteq\mathbb F_3^{13}\) are pairwise orthogonal;
if all three banks are occupied, each code is restricted to a
level of an affine quadratic whose homogeneous part is \(2Q\)
on its span, \(Q(v)=v\cdot v\).

For at most two banks, nondegeneracy of the ambient dot product
makes at least one span have dimension at most six. Hence the
combined word count is at most \(3^6+1449=2178<2187\).
With three banks, if all spans have dimension at most six the
target \(3\cdot3^6=2187\) is immediate. If one span has dimension
at least eight, the other two have dimensions at most five,
giving at most \(1449+2\cdot3^5=1935\).

It remains that \(U\) has dimension seven. The other two have
dimensions at most six. If both are at most five, the same
1935 bound applies. Otherwise, write one as \(V=U^\perp\)
of dimension six, and the other \(W\subseteq U\cap V\).
If \(Q|_V\ne0\), the \(V\) bank's quadric cuts every line in
a suitable parallel family at most twice, so it has at most
\(2\cdot3^5=486\) words. Since \(V\) is not isotropic,
\(\dim(U\cap V)\le5\), so \(W\) has at most 243 words.
The total is at most \(1449+486+243=2178\).

Finally suppose \(V\) is isotropic. Then \(V\subset U\) and
\(V\) is the radical of \(Q|_U\), with a one-dimensional
nondegenerate quotient. If the affine quadratic defining the
\(U\) bank has a nonzero linear part on \(V\), each parallel
line in one such direction has at most one selected word;
the \(U\) bank has at most \(3^6=729\) words. Otherwise the
\(U\) bank lies in at most two affine cosets of \(V\), from
a nonzero quadratic equation on the one-dimensional quotient.
If \(V\) contains a nonzero word of weight below six, that
weight is three because \(V\) is isotropic. Its parallel lines
allow at most \(3^5\) selected words per coset, so two cosets
have fewer than 729. If \(V\) has minimum weight at least six,
then two full distinct such cosets cannot be united at distance
six: they would give 1458 words, exceeding the published 1449
bound. Any conflicting pair across the two cosets supplies a
perfect matching of forbidden pairs by simultaneous translation
by \(V\). Thus at most 729 words can be selected across both.
Again the \(U\) bank's cap is 729, and with the other two,
the total is at most 2187.

## Smaller radii

The actual minimum squared separation \(\mu\) is at most 4
because of the within-period minimum. Below 4 the maximum
squared-distance levels to test are shown below; they follow
from the exact torus distance formula. For \(\mu\le2\) keep
at most three colors per word in each bank; above 2 keep
at most one. Ordinary ternary Hamming balls of radius 1 and 2
in length 13 have sizes 27 and 339 respectively, giving the
indicated code caps at distances 3–4 and 5. For distance 2 use
puncturing, and for 6 use the same published bound.

| \(\mu\) | Upper phase count |
| ---: | ---: |
| \(2/3\) | \(3^{15}\) |
| \(4/3\) | \(9\cdot3^{12}\) |
| \(2\) | \(9\cdot3^{10}\) |
| \(8/3\) | \(3\cdot3^{10}\) |
| \(26/9\) | \(3\lfloor3^{13}/339\rfloor\) |
| \(10/3\) | \(3\lfloor3^{13}/339\rfloor\) |
| \(32/9\) | \(3\cdot1449\) |

All rows give a smaller density than the unit-radius witness:
use the exact comparisons \(M(9\mu)^{13}<2187\cdot36^{13}\),
with the common covolume \(3^{15/2}\).
