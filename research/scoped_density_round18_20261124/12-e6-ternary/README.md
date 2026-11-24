# 12D ternary E6 × E6: additive motifs excluded for all independent scales

Consider the period \(E_6\times E_6\), with independent positive squared
block scales \((\alpha,\beta)\), and arbitrary motifs from its third-torsion
quotient. This proof excludes every **lattice** improvement over the
Coxeter–Todd center density \(1/27\) in this family. It reduces any
nonlinear improvement to a single equal-scale mask \((9/4,9/4)\) at
squared minimum 2. The mask still has a gap: its Fourier cap is 1215
versus the 923 sites needed. In particular no global bound is asserted.

## Closing the asymmetric vertex

At squared scales \((9/5,18/5)\), separation 2 is equivalent to
\(d(u)+2d(v)\ge10\) for nonzero differences, where \(u,v\) are the
horizontal and vertical differences in \(V=E_6/3E_6\) and \(d\) is
least squared lift norm. The radical \(R\) of the mod-3 bilinear form
has order three. The only residues of depth at least 10 are the two
nonzero radical elements, so horizontal sites at any given vertical
coordinate lie in one radical coset. Call a fibre *heavy* if it has
at least two sites, and light otherwise. A heavy fibre has at most
three sites.

For two heavy fibres all three horizontal difference lifts in their
quotient class occur, since two subsets of size at least two in a group
of order three have full sumset. Thus vertical depth 2 forces their
horizontal labels to differ by a nonzero isotropic element in the
nondegenerate quotient \(V/R\); vertical depth 4 forces unequal labels.
The least norm of a nonzero quotient difference is 2, 4, or 6
according to its quadratic value; this is checked by finite norm
enumeration.

Take the literal depth-2 six-clique \(S\) supplied in the replay.
If all six fibres there are heavy, the six distinct labels lie in an
affine totally isotropic plane \(P\). To see this, translate one label
to zero, use pairwise isotropy to obtain mutual orthogonality, and
use the Witt index bound two in the nondegenerate rank-five quotient.
Six points cannot lie on one affine line over \(\mathbb F_3\).

There are thirty additional sites \(U\), each with at least four
depth-2 neighbors in \(S\). A heavy site in \(U\) also has its label in
\(P\): four distinct labels span that affine plane, and the differences
of their isotropy constraints give orthogonality to its directions.
The remaining isotropy condition and maximality then put the label in
the plane. The thirty sites partition into three sets of ten, each
making a 15-clique of depths 2 or 4 with the same five seed sites. All
heavy members have distinct labels, at most nine in \(P\). Thus if
\(S\) is all heavy, each ten-set has at least six light sites and
\(U\) has at least eighteen.

Suppose the number of light sites is at most 127, and pad its set to
exactly 127 by declaring more sites light. The local implication above
still holds. Let \(X,Y\) be the light occupancies of a uniformly random
common translate and Weyl image of \(S,U\). The event \(X=0,Y\le17\)
is impossible by the local argument. The verifier supplies a rational
quadratic polynomial \(f(X,Y)\) at most the indicator of that event
at every integer pair in \([0,6]\times[0,30]\), with
\(\mathbb E f\ge1221965/238912254>0\), a contradiction.

Here the moments \(X^2,Y^2,XY\) are dot products of the orbit-averaged
light autocorrelation with the three corresponding tables of ordered
pair counts. For that autocorrelation \(p\) one has \(p\ge0\), \(p_0=q\)
and \(\sum_j |O_j|p_j=729q^2\), with \(q=127/729\). Its Fourier
transform is nonnegative, so \(Hp\ge0\) for the exact real Fourier
rows \(H\). The replay proves the expectation bound by a termwise
rational comparison using those constraints and the first moments
\(6q,30q\). No LP solver is needed.

Thus at least 128 fibres are light and the motif has at most
\(3\cdot729-2\cdot128=1931\) sites, below the 1935 needed. The period
covolume is \(3(9/5)^3(18/5)^3\), and the squared radius is \(1/2\).
This gives center density at most \(30171875/816293376<1/27\), the
[Coxeter–Todd density](https://doi.org/10.4153/CJM-1953-043-4).

## Full continuum of scales

Normalize the actual squared minimum to 2. The period forces
\(\alpha,\beta\ge1\). Clipping either scale to 9 preserves separation
since every nonzero least norm is at least 2. Swap blocks so
\(1\le\alpha\le\beta\le9\). For each fixed motif the admissible set
in this triangle is a polytope whose cutting lines are
\(a\alpha+b\beta=18\), with \(a,b\) least norms in the finite quotient.
The log of \(\alpha\beta\) is concave, so the minimum can be taken
at a vertex, and the replay checks all 64 possible arrangement vertices.

At a vertex the normalized, Weyl-averaged autocorrelation \(A\) in
product orbits is nonnegative, is 1 at the origin, and vanishes on
forbidden differences. Each row of \(H\otimes H\) dotted with \(A\)
is nonnegative by Fourier positivity: the product of real character
parts is the mean of two joint real character sums, with one character
kept or inverted. The exact rational dual at a vertex combines those
rows with nonnegative multipliers to bound the motif size. Beating
\(1/27\) requires more than \(64(\alpha\beta)^3/9\) sites. The duals
leave the asymmetric vertex closed above and the equal vertex
\((9/4,9/4)\), where the separation mask is \(d(u)+d(v)\ge8\)
and the cap is 1215. An additive motif has order a power of three,
so it has at most 729 sites there and cannot improve. A lattice that
is a union of these periodic cosets has an additive motif after
translation of one phase to zero.

Replay: python -B scales.py (Python standard library, exact only).

In the remaining equal-scale mask, the checked refinement
\(M\le30213/38+(2268/19)P_{\rm axis}/M\) counts ordered pairs
\(P_{\rm axis}\) differing only by a nonzero radical element in
one coordinate. It caps motifs without such pairs at 795, so an
improvement with \(M\ge923\) needs at least 990 ordered axis
radical pairs. See python -B equal_axis.py. This does **not**
reduce the unconditional equal-mask cap of 1215.
