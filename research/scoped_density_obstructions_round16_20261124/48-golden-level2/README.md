# Golden deformations of 2-modular lattices, and a level-10 theta bound

These are upper bounds within specified construction families, not
upper bounds on arbitrary lattice packings or new density records.

## Exact level-2/golden obstruction in rank 48

Let \(L\) be an even 2-modular lattice in dimension \(n\) divisible by
16. Write its original inner product as \(G\). Suppose \(S\) is
self-adjoint for \(G\), \(S^2=5I\), and \((I+S)/2\) preserves \(L\).
Allow **all** positive metrics \(G(aI+bS)\) for this action, not just
integral or modular twists.

The table gives rigorous upper bounds on their Hermite invariants.
There is no assertion that every upper value is attained.

| Rank | Upper Hermite invariant in this twist family |
| --- | --- |
| 16 | \(4/\sqrt2\) |
| 32 | \(6/\sqrt2\) |
| 48 | \(18/\sqrt{10}\) |
| 64 | \(10/\sqrt2\) |
| 80 | \(12/\sqrt2\) |
| 96 | \(32\sqrt{5/2}/5\) |

In particular the rank-48 family is strictly below the \(P_{48}\)
Hermite invariant of 6. The determinant of an even 2-modular rank-\(n\)
lattice in this convention is \(2^{n/2}\).

**Proof.** Put \(q=G(x,x)\), \(s=G(x,Sx)\), and
\(U=(I+S)/2\). The identity \(q(U^{\pm1}x)=(3q\pm s)/2\) implies
\(s\equiv q\bmod 4\) and \(|s|\le q\) on the shortest shell. Field-unit
reduction gives the fundamental range \(|t|\le1/5\) for \(I+tS\), and
its determinant-normalized squared lengths divided by the original
Hermite scale are \((q+ts)/\sqrt{1-5t^2}\). The real eigenspaces each
have dimension \(n/2\), since \(X^2-5\) is irreducible over the rationals.

The modular bound gives \(m=2+n/8\) as the upper bound on the original
minimum. If the parent has minimum at most \(m-2\), take an original
shortest vector and its unit translates. In logarithmic coordinates,
as in the rank-48 real quadratic analysis in Round 14, their minimum
is at most \((m-2)\cosh(\log((1+\sqrt5)/2))=(m-2)\sqrt5/2\)
in determinant-one coordinates. For clarity, for \(0\le r\le2R\),
\(R=\log((1+\sqrt5)/2)\), the two witness curves are
\(q\cosh r+z\sinh r\) and \(q\cosh(r-2R)+z\sinh(r-2R)\), \(z=s/\sqrt5\).
Their minimum is bounded by their intersection as functions of \(z\),
\(q\cosh R/\cosh(r-R)\), and the endpoints follow directly.
Divide by \(\sqrt2\) for the 2-modular starting determinant.

If the parent is extremal with minimum \(m\), its shortest shell is
a spherical 7-design. In this case one can see the required statement
directly: harmonic theta series of degrees two, four and six have
weight \(n/2+j\) on \(\Gamma_0(2)\), and both cusps have vanishing
order at least \(m/2\), since the scaled dual is isometric. There is
no assumption on the action of that isometry on a harmonic polynomial;
it merely changes the polynomial at the second cusp. The sum of the
two orders is at least \(m> (n/2+j)/4\) for \(j\le6\), contradicting
valence unless the series vanishes. Antipodality gives odd moments.
Consequently \(\mathbb E s=\mathbb E s^3=0\) and
\(\mathbb E s^2=5m^2/(n/2+1)\).

There are shortest-shell witnesses of both signs with absolute pairing
at least \(h=4,6,4,6,8,6\) in the table's six respective ranks.
For ranks 16 and 48 this follows from nonzero variance and the
grid of multiples of four. For ranks 32, 64 and 96, if the negative
side has no pairing \(\le-6\), the cubic \(s^3-4s\) is nonnegative on
the remaining grid and has zero mean, forcing variance at most 4;
use the opposite sign for the other end. In rank 80 the analogous
polynomial \(s^3-16s\) forces variance at most 16 unless there is a
witness at or below \(-8\). The actual variances exceed the respective
bounds. For each half of the \(t\)-interval use its witness to bound
the Hermite invariant by
\((m-h|t|)/\sqrt{2(1-5t^2)}\). Its derivative away from zero has sign
\(5m|t|-h\), so the maximum is at an endpoint. The table also covers
the weaker parent bounds above.

The obstruction includes the successive endpoints of any even
unimodular lattice with commuting self-adjoint integral actions
\(A^2=2I,\ S^2=5I\), with \(G(x,Ax)\) even and \((I+S)/2\) integral.
Indeed \(2I+A\) is an even 2-modular twist: its twice-dual Gram
operator is \(2I-A=(2I+A)(1+A)^{-2}\). The remaining action stays
self-adjoint and golden-integral in that metric. This is a substantial
class of mixed endpoints, not every strongly 10-modular lattice.

There is also a reverse test. Suppose an even strongly 10-modular
lattice has such an integral self-adjoint golden action \(S\), and
the factor-five similarity from \(L^\#\cap\tfrac15L\) onto \(L\)
is precisely \(S\). Suppose additionally its factor-two similarity
\(R\) commutes with \(S\). Here factor means the ratio of the squared
norms, and the duals refer to the original metric \(H\).
Then in rank 48 it cannot have minimum 20: it is a golden endpoint
of an even 2-modular lattice in the family above.

In fact, with \(C=(5I+S)/2=2I+U\), one has \(SC^{-1}=U^{-1}\).
Thus \(C^{-1}L=L^\#\cap\tfrac15L\) and \(G=HC^{-1}\) is integral
with determinant \(2^{n/2}\). It is even: for
\(\ell(x)=G(x,x)\bmod2\), evenness of \(H=G(2I+U)\) gives
\(G(x,Ux)\equiv0\bmod2\). Hence
\(\ell(Ux)=\ell(x)\), while \(U^2=U+I\) gives \(\ell(U^2x)=0\).
This forces \(\ell=0\). Also \(L^\#_G=CL^\#_H\).
The image of \(C(L^\#_H\cap\tfrac12L)\) modulo \(L\) has order
\(2^{n/2}\), because \(C\) is invertible modulo two (its determinant
is odd) and the original factor-two similarity gives that index.
It therefore generates \(L^\#_G\) together with \(L\).
Since \(R\) commutes with \(C\), it is a factor-two similarity
for \(G\) and maps these generators into \(L\); determinant
comparison gives equality. Consequently \(G\) is even 2-modular,
as needed.

## Separate bound for theta series invariant under full level 10

If an even rank-\(n\) level-10 lattice has theta invariant under all
Atkin–Lehner involutions, its theta lies in \(M_{n/2}(\Gamma_0(10)^+)\).
Use Gram determinant \(10^{n/2}\) when translating these norm bounds
to a strongly modular density. At ranks 64, 72, 80, 88 and 96 the
ordinary valence extremal minimum is impossible by an exact negative
next coefficient. The respective squared minima are at most
\(24,26,30,32,36\). At rank 48 the upper value is 20; a lattice attaining
20 must have 144,668,160 norm-20 and 157,248,000 norm-22 vectors.
This conditional form does **not** establish that such a lattice exists
or rule out every rank-48 level-10 candidate.

For a short exact calculation put \(e=(\sum_{d\mid10}d^2 E_4(dz))/130\)
and \(\phi=(10E_2(10z)+5E_2(5z)-2E_2(2z)-E_2(z))/12\).
Then \(f=13(\phi^2-e)/28\) has weight four and order one.
To see the invariance of \(\phi^2\), replace \(E_2\) by its modular
nonholomorphic completion. The common nonholomorphic terms cancel;
the normalized functions \(dE_2^*(dz)\) are permuted by the
squarefree-level involutions, and the signs in \(\phi\) form a
character of this four-element permutation action. Finally
\(g=\prod_{d\mid10}\eta(dz)^4\) has weight eight and order three.
It is the \(\Delta_{10}^3\) in Shigezumi's level-10 basis.
Alternatively, its eta-quotient cusp orders on \(\Gamma_0(10)\) are
all three; the involutions preserve its divisor, and the multipliers
are \(+1\) on the trace-zero representatives of \(W_5,W_{10}\), by
evaluation at their fixed points where the eta product is nonzero.
These two involutions generate the full extension.

For \(h_0=n/8\), monomials \(e^c f^a g^b\), \(c+a+2b=h_0\),
give all leading orders \(0,\ldots,\lfloor3n/16\rfloor\).
Valence on the full group bounds every holomorphic vanishing order
at infinity by \(3(n/2)/8\), so selecting one monic monomial
per order is a basis. Successively cancel all coefficients through
that largest order to obtain the unique putative extremal theta.
The next coefficients are computed exactly in verify.py. Its
negative coefficient at exponent \(3n/16+2\) rounded down for
noninteger \(3n/16\) excludes the extremal minimum in the listed ranks.
With one less prescribed zero, the same computation also gives exact
necessary lower counts on the first shell (see the replay).

## Reproduction and sources

Run python -B verify.py, using only the standard library. Definitions
of center and Hermite normalization agree with the [Cohn density table](https://cohn.mit.edu/sphere-packing/).
The two modular inputs are [Rains–Sloane, The Shadow Theory of Modular
and Unimodular Lattices](https://neilsloane.com/doc/Me231.pdf), and
[Shigezumi, On the zeros of certain modular functions for the normalizers
of congruence subgroups](https://citeseerx.ist.psu.edu/document?doi=0210841cffb47afe10ec0a971c1fe1e20719249f&repid=rep1&type=pdf),
section 10.1, formulas (196)–(198). Compare also [Nebe's account of
Venkov's theory](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf).
