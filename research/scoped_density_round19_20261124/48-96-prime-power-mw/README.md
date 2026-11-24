# The regulator cost defeats two prime-square Mordell–Weil constructions

Consider, over \(\overline{\mathbb F}_p(t)\), the sextic twist
\(y^2=x^3+t^q-t\) at \((p,q)=(5,25)\), and the quartic twist
\(y^2=x^3+(t^q-t)x\) at \((p,q)=(7,49)\).
Their geometric Mordell–Weil lattices have ranks 48 and 96.
In Shioda's height normalization the exact \(p\)-valuations of
their Gram determinants are 48 and 96, respectively.
Thus merely going to prime-square translation groups does not
produce competitive parents: their center density is at most
\(3^{-24}\) times the \(P_{48}\) value and \(3^{-48}\) times the
Bos–Conway–Sloane 96D value, respectively. These are bounds on
two named lattices, not on sections, neighbours, or all packings.
No catalogue lower bound changes.

The sextic case lies in the family of
[Griffon and Ulmer](https://arxiv.org/abs/1903.03901).
The point of the calculation below is a single integral
Frobenius-orbit argument that also handles the different quartic,
and distinguishes both from lower-cost prime-degree constructions.

## Tame characters and Frobenius divisibility

For \(k=6\) or 4 let \(C\) be the smooth curve \(v^k=t^q-t\), and let
\(E_0\) be the corresponding constant curve \(Y^2=X^3+1\) or
\(Y^2=X^3+X\). The twist becomes constant after \(x=v^2X,\ y=v^3Y\),
and its surface is birational to the diagonal quotient of \(C\times E_0\).
The diagonal action is \(v\mapsto\zeta v,\ (X,Y)\mapsto
(\zeta^{-2}X,\zeta^{-3}Y)\). Here the action on \(E_0\) written
with positive exponents has holomorphic character \(-1\bmod k\).

Let \(M=k(q-1)\); all its prime factors are prime to \(p\).
For a primitive \(M\)-th root of unity \(\xi\), the curve scaling
\((t,v)\mapsto(\xi^k t,\xi v)\) gives an integral eigenspace
decomposition of \(H^1_{\mathrm{cris}}(C)\) over a finite Witt ring
containing the tame roots. A holomorphic differential basis is
\[
 t^{i-1}dt/v^j,\quad
 1\le j\le k-1,\qquad 1\le i\le\lfloor(jq-1)/k\rfloor.
\]
Indeed at a simple finite branch point \(\operatorname{ord}_v(dt)=k-1\);
at the unique infinity, \(\operatorname{ord}(t)=-k\) and
\(\operatorname{ord}(v)=-q\), so the displayed range is exactly
the regularity condition there. Their number is the genus
\((k-1)(q-1)/2\). Their characters are \(ki-j\bmod M\).
The antiholomorphic quotient \(H^1(C,\mathcal O)\) has the negative
characters by duality. The combined list is exactly all \(a\bmod M\)
not divisible by \(k\), each once, which the exact replay checks.
Thus the crystalline character lines have rank one.

Write \(\epsilon(a)=1\) for the holomorphic characters and zero
for the others. Frobenius takes line \(a\) to line \(pa\bmod M\)
by multiplication by \(p^{\epsilon(a)}\) times a unit (with Witt
semilinearity). For clarity, in \(H^1\) the image of Frobenius
contains \(pH^1\) by Verschiebung, so the only possible valuations
are 0 and 1. Modulo \(p\) Frobenius kills the holomorphic
differentials, and its kernel has dimension the genus (the
Dieudonné \(H^1\) of the Jacobian). This gives the exact valuation
on each of the distinct tame lines.

In the diagonal-invariant cross term
\(H^1(C)\otimes H^1(E_0)\), the active \(a\) have
\(a\equiv 1,-1\bmod k\), matching the character of \(E_0\)
because its action on the quotient is inverse. The exponent of
\(F/p\) on the line \(a\) is therefore
\[
 d(a)=\epsilon(a)+\mathbf 1_{a\equiv-1\bmod k}-1 .
\]

## Integral fixed lattice and regulator

The replay checks that each Frobenius orbit of the *full*
\(H^1(C)\) has exactly half holomorphic positions.
Thus the Jacobian is supersingular, as is \(E_0\), and we may
pass to a common finite constant field over which their Frobenius
powers are the same scalar on \(H^1\). For instance use isogeny
to products of supersingular elliptic curves and take a further
even extension over which those isogenies are defined.
Also take the constant degree divisible by all the orbit lengths
and large enough that all geometric divisor classes are defined.

The integral comparison is important (a rational Tate space
alone would not determine a regulator). By
[Ulmer, Explicit Points on the Legendre Curve III, Theorem 5.2(1)](https://msp.org/ant/2014/8-10/ant-v8-n10-p04-s.pdf),
the \(p\)-completed cross Néron–Severi lattice of this product
identifies with the integral \(F=p\) crystal fixed lattice,
compatibly with the cup pairing. This is the \(p\)-adic
[Waterhouse–Milne integral Tate homomorphism theorem](https://www.jmilne.org/math/articles/1971a.pdf)
and Poincaré duality for the two factors. It commutes with
the tame diagonal action; taking its invariants is exact since
\(p\nmid k\). The product cross crystal has a perfect pairing,
and so does its invariant part (averaging is a self-adjoint
integral projection).

For one active orbit \(a_0,\ldots,a_{\ell-1}\), set
\(s_i=\sum_{j<i}d(a_j)\). Its total is zero because the slope is
one. The valuations allowed by integrality in the fixed lattice
are given by the minimal offsets \(s_i-\min_j s_j\).
To see that the offsets are attained without further index,
rescale the lines by those powers and solve the remaining unit
recurrence. The iterated unit cocycle has norm one because the
chosen finite-field Frobenius on the cross term is scalar.
Unramified Hilbert 90 gives a unit solution, and the first
coordinate then varies over the Witt ring of \(\mathbb F_{p^\ell}\).
A \(\mathbb Z_p\)-basis of this unramified extension has an
invertible matrix of conjugates after scalar extension
(the finite residue extension is separable). Thus the index
valuation after Witt scalar extension is exactly
\(\sum_i(s_i-\min_j s_j)\). Summing over active orbits gives
\(\sigma=24\) and 48. The corresponding cup discriminant has
valuation \(2\sigma\).

The same \(p\)-part holds on the elliptic surface. In the tame
diagonal quotient the pullback and descent on divisors differ
by multiplication by \(k\), a \(p\)-unit, and scale the
intersection form by \(k\). The fixed points on the product
are isolated; the exceptional chains resolving cyclic quotient
singularities have determinants dividing \(k\), also \(p\)-units.
Equivalently, over \(\mathbb Z_p\) one splits off these chains
and identifies the complement by pullback and pushforward of
divisor classes (local Weil denominators divide \(k\)).
The two fibre classes in the product have a unimodular pairing.
Birational passage between smooth models changes the Picard
lattice by unimodular blow-up classes. Shioda's determinant
formula now gives the same \(p\)-valuation for the MW Gram
determinant: both the fibre root determinants and the torsion
terms in these examples are prime to \(p\).
The MW rank equals the number of active lines, also obtained
directly by identifying equivariant maps from \(C\) to \(E_0\)
up to constants with the equivariant homomorphisms of Jacobians.

## Short sections and the density comparison

For \(p=5,q=25\) the Euler parameter is \(\chi=5\).
The 25 finite fibres of type II and the infinity fibre of type II*
have trivial component group, so a nonzero height is at least
10. The base change \(s=t^5+t\) satisfies \(s^5-s=t^{25}-t\);
the rational elliptic surface \(y^2=x^3+s^5-s\) has six type II
fibres and Mordell–Weil lattice \(E_8\). Pulling back a norm-two
vector gives a vector of norm exactly 10 (height scales by 5).
There is no nonzero torsion, since the fixed constants under
the order-six automorphism are trivial. The Gram determinant
is an integer divisible by \(5^{48}\), so the center density
is at most \((10/4)^{24}/5^{24}=2^{-24}\).

For \(p=7,q=49\) there are 49 type III fibres and one type III*
at infinity, making the product of fibre root determinants \(2^{50}\).
The torsion has order two, since a torsion section becomes a
constant fixed by the order-four automorphism after the tame base
change. Consequently the covolume is at least \(7^{48}/2^{24}\).

For a short vector take the degree-seven base change \(s=t^7+t\).
On \(y^2=x^3+(s^7-s)x\), set
\(A=s(s-1)(s-2)(s-3)\), \(B=(s-4)(s-5)(s-6)\).
Choose \(\lambda^2=-2\) over characteristic seven.
Then \(A+\lambda B\) is a square of a monic quadratic,
as checked exactly in the replay. Pick \(c^2=\lambda^{-1}\),
\(d^2=c^3\) in the algebraic closure, and use \(x=cA\),
\(y=dA\sqrt{A+\lambda B}\). This is a polynomial section
meeting the nonidentity component at just the four roots of
\(A\) and not at infinity (the leading quartic \(x\)-term is
nonzero there). It is disjoint from the zero section.
The squared height is \(4-4/2=2\) on this \(s\)-surface, so
its pullback has height 14. The center density in rank 96 is
at most \((14/4)^{48}2^{24}/7^{48}\), whose ratio to the BCS
center density \(3^{48}/2^{24}\) is \(3^{-48}\).

For the height, local components, and determinant convention see
[Shioda, Some Remarks on Elliptic Curves over Function Fields](https://www.numdam.org/article/AST_1992__209__99_0.pdf)
and [Shioda, Mordell–Weil Lattices and Sphere Packings](https://doi.org/10.2307/2374791).
Run python -B verify.py for the finite character arithmetic,
the literal short-section identity and the density comparisons.

The three remaining standard CM twists for the same choices of
\(p,q\) are treated in [OTHER_TWISTS.md](OTHER_TWISTS.md).

As a cross-check, Corollary 9.2(3) of Griffon–Ulmer establishes that
the sextic regulator is a power of \(p\) over sufficiently large
constants. Thus in that case the density \(2^{-24}\) is exact.
The use of Ulmer's earlier paper above is only Theorem 5.2(1);
it does not use the part (2) addressed by the
[published erratum](https://msp.org/ant/2014/8-10/ant-v8-n10-x04-Errata-PointsOnLegendreIII.pdf).
