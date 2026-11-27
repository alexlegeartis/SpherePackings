# Several planar sites per word also cannot escape the 27D gap

The [per-codeword phase result](PER_CODEWORD_PHASES.md) admits an
extension: assign every word *any finite number* of planar sites modulo
the common real period \(T\), provided all sites of that word use the
same whole horizontal high-digit parity coset. Each word can choose a
different parity, and it can be omitted entirely. The period still
must be a product with horizontal factor \(4D_{25}\). The same bound
\(8/(6+\pi)\) times the 27D Vardy incumbent holds.

Write \(W_b=8\sum_{w\in b}m_w/2^{14}\) for the normalized total number
of sites of bank \(b\), and \(W=\max W_b\). Here \(m_w\) counts the
normal sites assigned to word \(w\), and \(W\) is not assumed to be
at most eight. As before \(R\le4Wt^{27/2}/V\) with \(t=\mu/32\).
For any given word, all its normal sites and their period translates
are separated by at least \(\sqrt\mu\): one can hold the actual
horizontal point fixed. Disjoint planar disks therefore give
\(m_w\pi\mu/4\le V\) and thus \(W\le32V/(\pi\mu)\). In particular
\[
R\le (4/\pi)t^{25/2}.
\]
Suppose that \(R>c=8/(6+\pi)>7/8\). Since \(\pi>3\), it follows that
\(t^{25}>(21/32)^2>(24/25)^{25}\). This is enough for \(A_0=\mu-28>0\)
and \(A_0\ge4t^{27/2}\), because
\(378-(25/2)\cdot32\cdot(24/25)<0\) in the same derivative argument.

Take a certified constant-weight-eight three-dimensional subspace \(U\)
in the heavier bank. In each coset of \(U\), the same two-colour Voronoi
disk estimate applies to all \(m\) assigned planar sites at once,
including arbitrary multiplicity of a single label. Indeed there is
only one word of each label in that coset, so two different sites
of one label (including all period copies) have separation at least
\(\sqrt{\mu/A_0}\ge\sqrt8\) in normalized units. They cannot both cut a
positive-area cap on the same site's disk. Complementary labels still
cannot cut caps by the exact horizontal distances, regardless of the
number of copies. Each of the \(m\) sites thus has at most three opposite
cap neighbors. Writing \(g=\pi/8-1/4\), the area bound is
\(V\ge m A_0(1+g)\) in every coset. Averaging \(m\) over the cosets
gives \(V\ge W A_0(1+g)\). Now
\(R\le1/(1+g)=8/(6+\pi)\), contradiction.

The same literal subspace witnesses and checker apply. The additional
signs are the two elementary exact inequalities
\((21/32)^2>(24/25)^{25}\) and
\(378<(25/2)\cdot32\cdot(24/25)\).
