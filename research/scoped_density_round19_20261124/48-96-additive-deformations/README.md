# A regulator lower bound for linearized prime-square deformations

Let \(f=t^{p^2}+b t^p+c t\), with \(c\ne0\), over
\(\overline{\mathbb F}_p\). The four twisted Mordell–Weil parents
below are all uncompetitive, **for every** \(b,c\), not merely
for \(t^{p^2}-t\). The cubic 48D analogue is not closed here.

| Rank | \(p\) | Equation | Lower \(\frac12v_p(\det\mathrm{MW})\) | Center-density ratio cap |
|---:|---:|---|---:|---:|
| 48 | 5 | \(y^2=x^3+f^3\) | 18 | \(<0.929\) |
| 48 | 5 | \(y^2=x^3+f\) | 8 | \(<0.541\) |
| 96 | 7 | \(y^2=x^3+f^2x\) | 36 | \(<0.000003\) |
| 96 | 7 | \(y^2=x^3+fx\) | 22 | \(<0.118\) |

Comparators are \(P_{48}\) and BCS96. These are bounds for the
named geometric Mordell–Weil lattices, not for other carriers or
their sublattices. Full regulators of the deformed families are
not claimed. The exact constant-coefficient case is in
[the preceding regulator calculation](../48-96-prime-power-mw/README.md).

## Supersingularity and short sections persist

For the corresponding tame degree \(k\) put \(C:v^k=f(t)\).
Its translation group is the two-dimensional \(\mathbb F_p\)-space
of roots of \(f\). For any order-\(p\) subgroup \(H\), the quotient
is a curve \(v^k=s^p-\alpha s\) with \(\alpha\ne0\) up to scale:
take \(s=\prod_{h\in H}(t-h)\), and factor the linearized polynomial.
Over the algebraic closure rescaling gives \(v^k=s^p-s\).
The quotient map has degree \(p\). The Jacobian of \(C\) is
isogenous to the product of these \(p+1\) quotient Jacobians.
Indeed on \(H^1\) each nontrivial translation character occurs
in exactly the quotient for its kernel, and the invariant part
is zero since \(C/\ker(f)\) has coordinate \(v\) and genus zero.
Pullbacks of the quotient Jacobians thus induce a rational
cohomology isomorphism and an isogeny. Each is supersingular
by the prime-\(p\) character computation, so the integral
crystalline comparison of the preceding proof applies over
a sufficiently large finite field. The rank and short height
upper bounds are unchanged: use a degree-\(p\) pullback of the
same short section, with parameter and coordinates rescaled.
The simple zeros, degree, tame fibre types and torsion used
in the covolume comparisons are also unchanged.

## A uniform crystalline index lower bound from Cartier ranks

Keep the two elliptic Hodge lines as separate labels in the
diagonal-invariant cross term, including for \(k=2\), when their
characters restrict equally. The index defining \(\sigma\) is the
index of the integral \(F=p\) fixed lattice after Witt extension
inside the cross crystal. On any reduction modulo \(p\) thereof,
both conditions
\[
 F_\mathrm{cross}(z)\equiv0\bmod p,\qquad
 (F_\mathrm{cross}^2/p)(z)\equiv0\bmod p
\]
are necessary. The second is integral since the constant
supersingular elliptic curve over \(\mathbb F_p\) has \(F_E^2\)
equal to \(-p\) times coefficient Frobenius squared (trace zero).
The first condition on the elliptic antiholomorphic line has
codimension \(p^2-1-h_{k-1}\), where \(h_j=\lfloor(jp^2-1)/k\rfloor\);
the elliptic holomorphic line already has a factor \(p\).
On this latter line the second condition has codimension the
rank of \(F_C^2\bmod p\) in its matching character. On any
character the rank of \(F_C^2\bmod p\) is the rank of Hasse–Witt
on \(H^1(\mathcal O)\) there: \(F_C\bmod p\) kills the Hodge
filtration and injects its quotient into de Rham cohomology,
so composition has exactly that rank. Hasse–Witt is dual
to Cartier on the opposite character. These two conditions on
different source lines give independent codimensions, yielding
\(\sigma\ge p^2-1-h_{k-1}+\operatorname{rank}(\mathcal C|_{{j=k-1}})\).
Here Cartier interchanges \(j=1,k-1\), whose checked ranks
are equal. This is a lower bound on the index only.

## Finite Cartier rank formula

For \(t^{i-1}dt/v^j\), put \(j'=k-j\) and \(A=(pj'-j)/k<p\).
Cartier extracts from \(t^{i-1} f^A dt\) the terms of
degrees \(pn-1\), using \(v^{-pj'}\) to take out \(v^{-j'}\).
Write \(l,u,w\) for the counts of \(ct,t^{p^2},bt^p\)
in that power. There is only one \(l\) with \(0\le l\le A\)
and \(l\equiv-i\bmod p\). Set \(w=A-l-u\);
then the row is \(n=(i+l)/p+pu+w\), with multinomial coefficient
times \(c^l b^w\) before the Cartier \(p\)-th root.
For fixed \(i,n\) at most one term contributes. If \(b=0\),
replace \(c\) by 1 by column scaling. Otherwise take
\(\theta^{p-1}=b\). Then
\(b^w=\theta^{\,p(A-l)+(i+l)/p-n}\), a product of a row
and column scaling; \(c^l\) is also column-only.
The rank is therefore checked at \(b=0,1\) with \(c=1\).
The exact resulting ranks give the displayed bounds;
the cubic case gives only \(\sigma\ge14\), not enough
for its density comparison.

Run python -B verify.py. The local geometric inputs and integral
comparison are given in the preceding regulator packet; the
new replay checks the finite Cartier ranks and exact density ratios.
