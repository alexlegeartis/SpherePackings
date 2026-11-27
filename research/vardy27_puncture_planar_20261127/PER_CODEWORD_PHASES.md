# 27D: the gap survives arbitrary per-codeword planar phases

Keep the original [\(4D_{25}\) horizontal carrier and two stocks](README.md),
but now allow each retained word its own arbitrary planar site \(z_{b,w}\)
modulo a common full-rank real normal period \(T\). A site need not even
agree with those of other words with identical puncture label. This more
flexible family still has centre density at most
\[
\frac{8}{6+\pi}\cdot\frac1{\sqrt2}
\quad\text{(about 87.51% of the published Vardy incumbent).}
\]
Arbitrary deletion of whole codeword cosets and arbitrary high-digit parity
choices are permitted, as before. The horizontal/normal period remains a
product; neither retaining only part of a horizontal coset nor assigning a
single word multiple normal sites is allowed.

## The constant-weight simplex certificate

For each of the two banks and each triple \(F\) of puncture coordinates there
is a linear dimension-three subspace \(U\) of the stock whose seven nonzero
words all have weight eight and whose projection to \(F\) is an isomorphism.
See [literal witnesses](simplex_witnesses.txt) and
[integer checker](verify_simplex.py): every one of the 3276 triples is
verified for each stock independently. One convenient way to generate a
witness collection is to loop over weight-eight masks \(a<b<c\), taking
\(a\) smallest in its space, \(b\) smallest independent of \(a\), and \(c\)
smallest outside \(\langle a,b\rangle\). Accept when all seven nonzero XOR
sums have weight eight, and record the first accepted basis that covers
each full-rank column triple.

Every coset of \(U\) thus has one word of each label \(f\in\mathbb F_2^3\),
and every distinct pair in that coset has original distance eight. For the
retained words in any coset, the exactly realizable horizontal distances
are \(4(8-\mathrm{wt}(f+g))\) in squared units, by the sign/high-carry lemma
of the main proof. Hence these words have planar torus distances squared
at least \(A_0=\mu-28\) between opposite colours, at least \(2A_0\) within a
colour, and at least \(\mu-20>2A_0\) between complementary labels, when
\(A_0>0\). These implications compare only words in one coset, so there is
no longer any need to identify the phases of different words in a fiber.

## Averaging over cosets forces a dense geometric witness

Let \(t=\mu/32\), \(V=\operatorname{covol}(T)\), and
\(W_b=8|S_b|/2^{14}\), where \(S_b\) is the full retained set of bank \(b\).
As in the main proof, \(W=\max W_b\le8\) and the actual ratio \(R\) to the
27D incumbent is at most \(4Wt^{27/2}/V\). Let \(c=8/(6+\pi)>7/8\).
If \(R>c\), the lattice-area bound gives
\[
W>4c\sqrt3\,t^{-25/2}>6,\qquad
t^{25}>3c^2/4>147/256>(39/40)^{25}.
\]
So \(A_0>0\), and \(A_0\ge4t^{27/2}\) by the same monotonicity as in the
main proof, since \(378-(25/2)\cdot32\cdot(39/40)<0\).

In the heavier bank the mean retained occupation of a coset of \(U\) is
exactly \(W\). Thus a coset with at least seven retained words exists.
Use any seven of those sites on the torus. In units where lengths are
divided by \(\sqrt{A_0}\), the main proof's Voronoi disk argument bounds
the area at each site by one from below, with an additional
\(g=\pi/8-1/4>1/8\) for every one of the four potential caps that is
missing. Only opposite-colour neighbors at distance below \(\sqrt2\)
can create caps. These neighbors must have distinct labels, since the
period has normalized minimum squared length at least \(\mu/A_0\ge8\).
The \(4+3\) colour counts among seven sites already yield four missing
caps. Exactly three complementary pairs survive among seven labels; none
can produce caps at either endpoint, yielding six more missing caps.
Consequently
\[
 V\ge (7+10g)A_0,\qquad
 W>c(7+10g)>(7/8)(7+10/8)>7.
\]
Average again to get a *full* coset of eight retained labels. Each of its
eight sites lacks at least one cap, due to its complementary label. Now
\(V\ge8(1+g)A_0\), and the ratio bound gives
\(R\le1/(1+g)=8/(6+\pi)=c\), contradiction.

The cap estimate used here is pure geometry, and the simplex witnesses
are sufficient regardless of the cross-bank constraints. Both banks are
tested because the heavier one is not known in advance. The result is
not an upper bound outside the stated construction model.
