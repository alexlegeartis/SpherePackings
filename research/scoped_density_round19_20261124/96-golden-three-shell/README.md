# 96D golden branch: both marked Jacobi parities couple norms 10 and 12

Suppose an even unimodular rank-96 lattice has squared minimum 10
and an integral self-adjoint endomorphism \(U^2=U+I\). Write
\(S=2U-I\), \(s(y)=y\cdot Sy\). This is the \(A=0\) branch of
[round 18](../../scoped_density_round18_20261124/96-golden-anchored/README.md).
Its known signed shell has \(1,655,821,440\) shortest vectors of
\(s(x)=-10\), so choose any one. Then \(x,Ux\) are orthogonal
norm-10 anchors. This note uses exactly this same pair throughout.

## The exact new exclusion of the old witness

Consider
\[
 \Theta(\tau,z,w)=\sum_{y\in L}
 q^{y^2/2}r^{x\cdot y}t^{Ux\cdot y}.
\]
It has weight 48 and diagonal Jacobi index \((5,5)\). Central
inversion leaves only its even-even and odd-odd parity parts.
By [Wang–Williams, Theorem 2.4](https://arxiv.org/abs/2102.11190),
for the even-even part use the tensor products
\(\phi_{-2,1}^a\phi_{0,1}^{5-a}\) in each variable with elliptic
coefficients. Their leading Laurent polynomials are powers of
\(r+r^{-1}-2\) and \(r+r^{-1}+10\); in the coded coefficient basis
the \(36\) bivariate products form an invertible square matrix.
Successively match \(q^0,\ldots,q^4\) to \(1,0,\ldots,0\), at each
step using \(\Delta^j E_4^uE_6^v\) with \(v=0\) or 1 to obtain
the needed elliptic weight (a missing weight-2 polynomial is
never requested by a nonzero coefficient). The remainder is
divisible by \(\Delta^5\); in the quotient the elliptic weight
for index exponents \(a,b\) is \(2(a+b)-12\). Only \(a+b=6,8,9,10\)
remain, eleven products in total. Their \(q^5\) coefficient arrays
are independent, so all subsequent coefficients are determined by
the marked norm-10 array. The odd-odd remainder uses in each
variable an extra \(\phi_{-1,2}\) with leading \(r-r^{-1}\),
and powers \(\phi_{-2,1}^a\phi_{0,1}^{3-a}\); the quotient
elliptic weights \(2(a+b)-10\) allow only \((2,3),(3,2)\).
Again their two distinct leading Laurent products determine
later coefficients.

The exact round 18 rational witness does satisfy the **full** new
even-even \(q^5\) identities. But the uniquely forced \(q^6\)
coefficient is negative at absolute dot labels \((4,5),(5,4)\)
(equal to \(-13781219850/7471\)) and \((5,5)\) (equal to
\(-50323730375/7471\)). This is impossible for the even-even
symmetrization of an actual nonnegative count. Thus this
specifically eliminates that formal witness, not every lattice.

## A new feasible exact cross-shell relaxation

The new model keeps one rational count for each orbit of
\((p,c,s)=(x\cdot y,Ux\cdot y,s(y))\) at each of squared norms
10 and 12 (222 and 441 orbits). In each orbit the Gram of
\((x,Ux,y,Uy)\) is PSD and the checked integer combinations
have no positive squared norm below 10. For each signed
shell use its known population and split covariance from
the Hilbert-theta identities of the previous proofs.
The ordinary shells are 11-designs by harmonic elliptic theta
series. In addition to each shell's negation, bridge maps
of the **sets counted** are literal bijections: include
\(y\mapsto y+ux+vUx\) for \(u,v\in\{-1,0,1\}\) when it has
the new specified norm, and \(y\mapsto Uy\) or \(U^{-1}y\)
when it changes norm 10 to 12. These are count equalities
between individual labels, not angularly averaged shell totals.
Both Jacobi parities then link the ordinary marked array at
norm 10 to the ordinary marked array at norm 12 coefficient by
coefficient. The replay also keeps the previous norm-10
harmonic-degree-14 rows and the four exceptional anchor vectors.

There is a **formal nonnegative rational witness** for all 652
checked rows, supported on 151 orbits including the exceptional
anchor orbit. Every positively used orbit additionally passes
a complete rank-four integer span check through squared norm
16, including exclusion of signed shell values whose population
is zero at \(A=0\). This checks valid support for this witness,
not an existence theorem for a lattice. In particular feasibility
with the same marked pair, these two shells, and all their Jacobi
cross constraints does not certify feasible joint choices
for unrelated anchors.

Replay: python -B verify.py (Python 3 + SymPy for small exact
linear algebra; no optimizer). Old replay files and old witness
are copied intact as standalone exact inputs.

## All affine-unit relations between the same two shells

The elementary 652-row witness has yet another exact failure.
The map \(y\mapsto U^{-1}y-x\) bijects norm-10 label
\((-4,2,2)\) with norm-12 label \((-4,-4,-4)\);
the two counts in the witness are 0 and 227,401,200.
This is an exact necessary condition from affine composition,
not just from counting either single step on its own.

The replay adds **all** maps \(y\mapsto U^k y+ux+vUx\)
that take classes in these two shells into either shell.
For a vector outside the anchor plane, the squared length of
its projection perpendicular to that plane grows in both
tails of the unit orbit (both golden eigenspaces must occur,
as the projection has rational lattice coordinates and the
eigenvalues are irrational). For projected squared norm
\(\rho\) and signed quadratic value \(\sigma\), the unit
recurrences are
\[
 (\rho,\sigma)\longmapsto
 ((3\rho\pm\sigma)/2,(\pm5\rho+3\sigma)/2).
\]
Thus a tail can stop permanently once \(\rho>12\) and its
outward-signed \(\sigma\) is nonnegative. For each unit
translate retained, list every translation in the anchor
plane whose two new dot magnitudes are at most 10; larger
magnitudes cannot give squared norm at most 12. Inside the
anchor plane, norm 10 has only the four exceptional labels
already fixed in one orbit, and norm 12 would require a sum
of two integer squares equal to 120, which is impossible.
There are 221 distinct additional pair identities or zeros
(including repeated restrictions already detected by prior rows).
The final exact 153-positive-orbit witness passes all **873**
rows and the full small-span arithmetic check. No integer
or full-lattice feasibility claim is made.

## Norm 14: a third shell, still jointly feasible in the relaxation

Both parity components at norm 14 are likewise predicted from
the single marked norm-10 array. The replay now compares those
coefficients to a third signed shell (728 orbits) and includes
every elementary shift from norms 10 or 12 to 14 and golden
unit move between those shells. It uses exact norm-14 signed
populations, split covariances and ordinary moments, and checks
the positively used norm-14 orbits exhaustively through squared
norm 16 in their four-vector spans.

The same exact 153-orbit affine-refined two-shell witness
extends, **without changing any of its counts**, by 145
additional positive norm-14 orbits. All 1953 three-shell
rows pass, as do the 221 complete two-shell affine rows.
We did not impose every long affine-unit composition involving
norm 14, nor any joint distribution across unrelated anchors.
Near-infeasibility or even infeasibility statuses from numerical
LP were unreliable under poor scaling in this calculation;
only the exact rational check is a result.
