# 96D golden branch: the single-marked-pair relaxation survives four shells

This is a checked **non-exclusion** within the \(A=0\) golden branch,
not a proposed lattice. Keep the same hypothetical extremal even
unimodular lattice with squared minimum 10, golden self-adjoint
operator \(U^2=U+I\), and orthogonal norm-10 anchors \(x,Ux\)
as in [Round 19](../../scoped_density_round19_20261124/96-golden-three-shell).
The previous 298 rational orbit weights extend **without changing
any value** to squared norm 16, with 258 additional positive orbit
weights. In particular a search restricted to the following verified
single-pair model cannot disprove this branch.

The first three marked signed shell arrays, including both parities of
the rank-two Jacobi form, remain fixed. The two parity decompositions of
[Wang–Williams, Theorem 2.4](https://arxiv.org/abs/2102.11190) give the
norm-16 array from the same norm-10 input by expanding the identical
weight-48, index-(5,5) Jacobi polynomial spaces one coefficient further.
The code checks all produced dot labels. At norm 16 use the [previous
golden Hilbert-theta signed populations](../../scoped_density_round17_20261124/96-product-two-point)
and their separate split covariance, together with the ordinary
rotational moments through total degree 10. Their total independently
matches the unmarked elliptic theta coefficient.

All positively used marked four-vector spans are checked up to squared
norm 16 with exact inverse-Gram coefficient bounds. In the new shell an
outside-anchor-plane vector gives full rank: the stable rational span
of \(x,Ux,y,Uy\) is a vector space over \(\mathbb Q(\sqrt 5)\).
Within the anchor plane the lattice is exactly \(\mathbb Zx+\mathbb ZUx\):
subtracting a nearest point of this orthogonal lattice would otherwise
produce a nonzero lattice vector of squared norm at most 5. So the
only plane vectors in the four checked shells are the four anchors
themselves.

The checker also adds **all** affine-unit bijections
\(y\mapsto U^k y+a x+b Ux\) between marked strata in any of the four
shells. The perpendicular squared length tends to infinity in both
directions of the golden orbit, giving an exact stopping check on each
tail. All possible horizontal translations with both dot magnitudes
at most 10 are retained, as follows from separation at distance
\(\sqrt {10}\); all four anchor exceptions are explicitly included.
The 3,257 affine identities therefore cover compositions, not merely
unit and translation moves taken separately.

Run `python -B verify.py` (Python and SymPy only) in this directory.
The optional argument is the earlier three-shell packet's path.
It verifies 2,166 equality rows, 3,257 additional affine checks,
and all support arithmetic without a floating-point optimizer.
A useful next constraint should link **different** marked anchor
pairs or exploit global integrality, rather than scaling the same
single-pair model indiscriminately.
