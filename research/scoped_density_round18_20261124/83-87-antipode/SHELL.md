# The norm-8 shell survives all small-determinant minimum normal triples

Let \(L\) be any **even** unimodular lattice of rank 88 and squared
minimum 8. For any independent norm-8 normal triple with Gram
determinant at most 384, at least **1,033,559** minimum vectors
(including both signs) are perpendicular to the triple. Thus the
85-dimensional period in this subspace has squared minimum exactly 8.
No assumption of primitivity of the displayed triple is needed.

The shell size is 168,498,000 by the extremal theta series
([Sakamoto, Proposition 6.2](https://arxiv.org/html/2609.15633v2#S6)).
Weighted theta in harmonic degree \(d=2\) or \(6\), divided by
\(\Delta^4\), has weight \(44+d-48\), namely \(-2\) or \(2\).
It vanishes, so both harmonic means on that shell vanish.

Write \(G=(v_i\cdot v_j)\) and \(b_i=v_i\cdot x\) for a shell point.
Apart from the six signed \(v_i\), each \(b_i\) is an integer in
\([-4,4]\) by minimality applied to \(x\pm v_i\). Also
\(b^\mathsf{T}G^{-1}b\le8\). The exceptional labels are the signed
rows of \(G\). Signing and permuting reduces \(G\) to sorted
absolute off-diagonals all positive or all negative (one choice
suffices with zeros). In the all-negative case their absolute sum
is at most 8 by minimality of \(v_1+v_2+v_3\). This leaves 26
candidate Gram types up to the determinant limit.

On polynomials in the labels the Euclidean Laplacian is
\(\sum_{i,j}G_{ij}\partial_i\partial_j\).
For a homogeneous polynomial \(P\) of degree \(d\), its harmonic
projection is \(\sum_j c_j\|x\|^{2j}\Delta^jP\) with \(c_0=1\) and
\(c_j=-c_{j-1}/[2j(88+2d-2j-2)]\). This follows by applying the
Laplacian to the expression and canceling successive terms.
Substitute \(\|x\|^2=8\) on the shell and form a rational combination
\(p_G\) of the resulting monomial projections of degrees 2 and 6,
plus a constant \(p_0\). The replay checks \(p_G(b)\le
\mathbf 1_{\{b=0\}}\) on all allowed ordinary labels.
By the harmonic means, the number of zero labels is at least
\(168498000p_0-2\sum_i p_G(G_{i,:})\). Its minimum over the 26
charts is strictly above 1,033,558. The finite checks do not use
the numerical optimizer that proposed the rational coefficients.

For a primitive normal lattice \(K\), projected phases depend
only on their class in \(K^*/K\). A literal norm-8 vector must
use a nearest normal representative in its class, since otherwise
subtracting \(K\) would contradict minimality. Thus for class \(c\)
of least squared normal norm \(q_c\), such a witness gives
horizontal distance squared \(8-q_c\); if no norm-8 vector has
the class, the distance is at least \(10-q_c\). The latter absence
is not proved by this shell-at-zero certificate.

Replay: python -B verify_shell.py (standard library, exact).
