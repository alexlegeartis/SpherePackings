# Strict local translation rigidity for the two near-record norm-12 sheets

Use the [47D primitive norm-12 period](README.md) and the two adjacent
canonical phases. At the fixed period, every sufficiently small nonzero
relative translation of the two sheets **decreases** their squared
minimum \(71/12\). This does not address collective changes of the
period metric or large noncanonical shifts.

Let \(A=\{x\in\Lambda:\|x\|^2=6,\ x\cdot v=1\}\); it is nonempty by
the moment calculation in the main note. Its horizontal projections
are exactly the contact vectors for one direction between those
two phases. Put \(a=x\cdot v\). On the whole norm-six shell
\(|a|\le6\). For each horizontal vector \(z\), apply the 11-design
identity to the degree-10 polynomial
\[
  (x\cdot z)\,a(a^2-4)(a^2-9)(a^2-16)(a^2-25).
\]
Its spherical average is zero by reflection of the horizontal
coordinate along \(z\). The terms at \(a=\pm2,\ldots,\pm5,0\) are zero.
At \(a=6\), \(x\mapsto v-x\) is an involution of the same bin and
negates the horizontal component, so this bin also has zero sum.
The \(a=-6\) bin behaves likewise. The two terms from \(a=\pm1\)
contribute twice the same multiple of the \(a=1\) horizontal sum.
It follows that the projected contacts have centroid zero.

They span the entire perpendicular space. For a unit horizontal
\(z\), use the degree-10 polynomial
\[
 F(a^2)(x\cdot z)^2,\qquad F(u)=-u(u-4)(u-9)(u-16).
\]
Among integral labels \(|a|\le6\), the only positive \(F(a^2)\)
occurs at \(|a|=1\). The spherical average of this product is
\(6804/325>0\); it follows by applying the coefficients
\((0,576,-244,29,-1)\) to the exact spherical moments
\[
 M_0=1/8,\qquad M_{j+1}=M_j\,\frac{72(2j+1)}{50+2j},
 \quad M_j=\mathbb E[a^{2j}(x\cdot z)^2].
\]
The design equates this average to the shell average. Thus some label-one
projected contact has nonzero component in every nonzero direction.

These finite contact vectors span and have zero centroid, so zero
lies in the interior of their convex hull: a supporting hyperplane
through zero would force every vector onto that hyperplane. Hence
there is \(c>0\) such that for every nonzero relative motion \(u\),
one contact \(y\) satisfies \(y\cdot u\le-c\|u\|\). Its new squared
length is at most \(71/12-2c\|u\|+\|u\|^2<71/12\) if
\(0<\|u\|<2c\).

For a separate useful metric diagnostic, the contact outer products
span the full symmetric matrix space as well. If \(T\ne0\) is a
horizontal real symmetric form, set \(G(x)=(x^{\mathsf T}Tx)^2\) and
use \(a^2(4-a^2)(36-a^2)G(x)\). This degree-10 polynomial has positive
factor in \(a\) only for \(|a|=1\). For
\(J_j=\mathbb E[a^{2j}G(x)]\), orthogonal Gaussian moment separation
and division by the spherical radial moment give
\(J_{j+1}/J_j=72(2j+1)/(52+2j)\). Thus
\(J_2/J_1=4,\ J_3/J_2=45/7\), and the spherical average is
\(144J_1-40J_2+J_3=(68/7)J_1>0\). This proves that no nonzero
\(T\) vanishes on all contact outer products. Spanning is not
positive eutaxy and by itself does not prove metric rigidity.
