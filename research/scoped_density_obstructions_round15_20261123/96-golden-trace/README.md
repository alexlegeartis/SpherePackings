# Golden trace deformations: exact 96D split and 5-modular reverse test

This is an exact family result, not an improved packing bound. No
96-dimensional parent in the high-density branch is constructed here.

Let \(L\) be even unimodular with inner product \(G\), and \(S\) be
self-adjoint, \(S^2=5I\), with \(U=(I+S)/2\in\operatorname{End}(L)\).
Thus \(U^2=U+I\) and \(U\) is an integral automorphism. Write
\(q=G(x,x)\), \(s=G(x,Sx)\). After scale and pullback by a power of \(U\),
every positive trace deformation \(aI+bS\) is represented by \(I+tS\)
with \(|t|\le1/5\): in the determinant-one form
\(\cosh(r)I+\sinh(r)S/\sqrt5\), pullback shifts \(r\) by
\(2\log((1+\sqrt5)/2)\). The Hermite invariant is the minimum of
\(q+ts\) divided by \(\sqrt{1-5t^2}\).

## Complete trace curve

If the parent has rank 96 and minimum 8, the **exact** unnormalized
minimum on the fundamental interval is \(8-h|t|\):

| Pairings on the shortest shell | \(h\) | Maximum Hermite invariant | Minimum in \((5I+S)/2\) metric |
| --- | ---: | ---: | ---: |
| All vanish | 0 | \(4\sqrt5\) | 20 |
| Nonzero, but no absolute value 8 | 4 | \(18/\sqrt5\) | 18 |
| Absolute value 8 occurs | 8 | 8 | 16 |

For extremal even unimodular parents in ranks \(24,48,72,96\) (minimum
\(m=4,6,8,10\), respectively), the curve is \(m(1-|t|)\). In rank 96
this is conditional on the existence of such an extremal parent.
Their 5-modular endpoint minima are \(2m\).

The unit gives
\[
(q(Ux),s(Ux))=((3q+s)/2,(5q+3s)/2).
\]
The inverse-unit norm is \((3q-s)/2\). As both norms are even and at
least the parent minimum \(m\), \(s\equiv q\bmod4\) and \(|s|\le3q-2m\).
On the shortest shell \(|s|\le m\); if \(s=-m\) occurs at \(q=m\), applying
the unit gives \(s=m\) at the same norm, and applying the inverse does
the opposite. For \(n=96,m=8\), the shell has mean \(s=0\): the
degree-two weighted theta series for \(s\) is divisible by \(\Delta^4\)
with residual weight two. The polynomial \(s\) is harmonic since
the two eigenspaces have equal dimension (the characteristic polynomial
of the rational matrix \(S\) is a power of \(X^2-5\)).
This gives the listed extreme pairings \(\pm h\).
For \(h=8\), \(q+s/5\ge(2q+16)/5\ge32/5\). For \(h=4\) the same
bound at \(q\ge10\) gives \(q+s/5\ge36/5\). For \(h=0\) it gives
\(q+s/5\ge8\) at \(q\ge12\); at \(q=10\), the only violating grid
value is \(-14\), which maps to \((q,s)=(8,4)\), impossible.
Use \(-S,U^{-1}\) for the opposite sign. Affine interpolation with
\(t=0\), and the shortest-shell witnesses, give exact equality.
The derivative of the normalized curve has sign \(40t-h\) on the
positive half, so the maximum is one of its two endpoints.

For the extremal parents, existence of \(|s|=m\) suffices: the unit
gives both signs and \(|s|\le3q-2m\) makes every other norm no shorter
throughout the interval. In ranks 24/48/72 the shortest shell is an
11-design, by the usual weighted theta argument. Projection onto one
of the two eigenspaces gives moments \(5m^2/(n/2+1)\) and
\(75m^4/((n/2+1)(n/2+3))\) of degrees two and four in \(s\).
The allowed grids force \(|s|=m\): respectively the variance is
positive, the variance exceeds 4, or the fourth moment exceeds
16 times the second. Rank 96 with minimum 10 follows from the
norm-10 calculation below at \(A=0\).

For nonextremal parents in ranks 24/48/72, and parents with minimum at
most six in rank 96, the endpoint minimum is smaller than the respective
ceiling stated above: the two endpoint forms are isometric by \(U\),
so one can take the smaller of a shortest vector's two endpoint norms,
giving at most \(5m/2\), rounded down to even.

## Equivalent search on the 5-modular side

Start instead with any even integral lattice \((M,H)\) and an
\(H\)-self-adjoint endomorphism \(S\) with \(S^2=5I\),
\((I+S)M\subseteq2M\) and \(SM^\#=M\), where \(M^\#\) is the \(H\)-dual.
The last condition specifies the self-adjoint 5-modular similarity;
an arbitrary 5-modular similarity is not assumed to have this property.
In rank 96 this class has minimum at most 20; equality holds iff the
reverse parent has minimum 10, or minimum 8 and full shortest-shell
cancellation.

Indeed, set \(C=(5I+S)/2=2I+U\). It is positive and \(SC^{-1}=U^{-1}\).
From an even unimodular parent, \(H=GC\) is even since
\(G(x,Ux)=G(Ux,Ux)-G(x,x)\) is even; the new dual is \(C^{-1}L\), which
\(S\) takes onto \(L\). Conversely \(G=HC^{-1}\) is integral
unimodular since \(C^{-1}M=S^{-1}U^{-1}M=M^\#\).
For evenness let \(\ell(x)=G(x,x)\bmod2\), a linear functional.
Evenness of \(H\) gives \(G(x,Ux)=0\bmod2\), hence
\(\ell(Ux)=G(x,(U+I)x)\bmod2=\ell(x)\).
But \(\ell(U^2x)=\ell(Ux+x)=0\), so \(\ell=0\). Apply the rank-96
even unimodular minimum bound of 10 and the exact curves above.

## Necessary shell counts in the high-density branch

Suppose there is no norm below 8 and all norm-8 vectors have \(s=0\).
Let \(A\) count these vectors (possibly zero for minimum 10).
Norm 10 allows only \(s=\pm2,\pm6,\pm10\) by the unit inequalities
and exclusion of \(\pm14\) as above. The ordinary weight-48 theta
series is the formal minimum-10 theta series plus \(A\Delta^4\).
Thus the total norm-10 count is \(565866362880-96A\).

Let \(H_k\) be the homogeneous harmonic in \(q,s\) led by \(s^k\),
with other monomials \(q^{2i}s^{k-2i}\), of degree \(2k\) in \(x\).
Its coefficients follow by canceling successive Laplacian terms:
\[
\Delta_{\mathbb R^{96}}(q^a s^b)
=2a(2a+94+4b)q^{a-1}s^b+20b(b-1)q^{a+1}s^{b-2}.
\]
In particular \(H_2=s^2-(5/49)q^2\) and
\(H_4=s^4-(30/53)q^2s^2+(25/901)q^4\).
The weighted theta series are multiples of
\(\Delta^4E_4,\Delta^4E_4^2\), so their norm-10 coefficients are
144 and 384 times their norm-8 coefficients. These identities and the
total give three independent equations, yielding

| \(\lvert s\rvert\) | Norm-10 count |
| ---: | ---: |
| 2 | \(462780864000-40A\) |
| 6 | \(99773856000-60A\) |
| 10 | \(3311642880+4A\) |

The odd-in-\(s\) weights \(H_1,H_3,H_5\) give zero at norm 10:
the series starts after order four, with residual weights \(2,6,10\)
after dividing by \(\Delta^4\), and none has a cusp form. The resulting
first, third and fifth moments vanish, giving equality of the positive
and negative counts within each row by a Vandermonde. This does not
require a Galois isometry. Nonnegativity gives \(A\le1\,662\,897\,600\).

The endpoint norm-20 shell consists exactly of
\((q,s)=(8,0),(10,-10),(12,-20)\) by the unit inequality. The last type
is in bijection with the first by \(U\). Thus the endpoint's kissing
number is exactly **\(1\,655\,821\,440+4A\)**. Its center density
\(5^{24}\) would be about 5.9 times the repository lattice comparator,
but these necessary conditions do not establish that a parent exists.

## Reproduction and sources

Run python -B verify.py (Python standard library only) to replay the
shell counts, harmonic polynomials, sign equations and endpoint
arithmetic. The proof uses the classical ring \(\mathbb C[E_4,E_6]\),
modularity of harmonic theta series for even unimodular lattices, and
the even unimodular minimum bound.

- [Previous trace-twist analysis](../../scoped_density_obstructions_round14_20261123/48-72-trace-twists/README.md)
- [Rains–Sloane, The Shadow Theory of Modular and Unimodular Lattices](https://neilsloane.com/doc/Me231.pdf)
- [Dursthoff, Extremal Lattices and Hilbert Modular Forms](https://publications.rwth-aachen.de/record/685470/files/685470.pdf), Chapters 5–6

## Transfer to 88D: a distinct non-Nebe test

The same exact \(8-h|t|\) curve holds for every rank-88 even unimodular
minimum-8 parent with the full golden action: the only added input was
mean-zero \(s\), and the degree-two theta series now vanishes already
by weight \(46<48\). If \(h=0\), the 5-modular endpoint has squared
minimum 20. Its center density \(5^{22}\) is about **36.2 times** the
repository's rank-88 lattice comparator and about **27.7 times** its
periodic comparator. In the other two branches it cannot improve the
current rank-88 lattice comparator.

The ordinary rank-88 series has weight 44 and four basis elements
\(\Delta^j E_4^{11-3j}\) for \(0\le j\le3\). Its norm-8 and norm-10
coefficients are 168498000 and 2480127344640. If \(h=0\), the weighted
degree-four and degree-eight series are multiples of \(\Delta^4\)
and \(\Delta^4E_4\), giving factors \(-96,144\) in place of the
96D factors. The same calculation gives the norm-10 counts
\(1970971200000,486367200000,22788944640\) at \(|s|=2,6,10\).
The odd-in-\(s\) weights again prove equality of the signs, since
the residual weights are \(-2,2,6\). The endpoint kissing count is
then \(11731468320\).

Gauvin's cyclotomic ideal lattices
[X4049](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/dim88e.html)
and [X4969](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/dim88f.html)
give concrete parents to test. Each is reported as even unimodular
with minimum 8. The ring of \(\mathbb Q(\zeta_{115})\) contains the
golden integers acting on the ideal; multiplication by the real
\(\sqrt5\) is self-adjoint for its trace inner product. No claim is
made here that either shortest shell fully cancels.

## Next-shell constraints in 96D

The same method determines the full absolute-value split at norm 12,
still under shortest-shell cancellation:

| \(\lvert s\rvert\) | Norm-12 count |
| ---: | ---: |
| 0 | \(1218201435360000+3212A\) |
| 4 | \(1453085210304000-1540A\) |
| 8 | \(286338229632000+2480A\) |
| 12 | \(14383631808000+370A\) |
| 16 | \(99773856000-60A\) |
| 20 | \(2A\) |

The unit bijects the final type to twice the norm-8 shell (one sign for
each of \(U,U^{-1}\)), and the penultimate to the norm-10 absolute
pairing-six type. For the four unknown counts use ordinary theta and
\(H_2,H_4,H_6\). The \(H_6\) theta has basis
\(\Delta^4E_4^3,\Delta^5\), and the first two coefficients were fixed
by the earlier shell split. The equations have nonsingular weighted
Vandermonde matrix, and the verifier checks their constant and linear
parts separately. The large \(s=0\) count gives an additional target
for geometric shell bounds, not a nonexistence theorem.
