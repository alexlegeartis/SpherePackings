# Arbitrary positive phase counts at normal diameter squared below 3/4

This strengthens the [six-label certificate](README.md) without changing
the rank-five Ne64 normal family. No finite cluster in the projected dual
of any of its admissible graphs beats the current 59D record **on the
ambient-minimum-minus-diameter certificate** when its squared diameter
is below \(3/4\). This still does not compute the true horizontal minima.

Write \(L_g\) for the full six-by-six weighted graph Laplacian. A dual
label difference is an arbitrary integer vector \(v\) on its six vertices
with zero sum, its five unanchored entries giving the original label.
Its squared normal length is \(3v^{\mathsf T}L_g^+v\). Since every degree
is 12, the Laplacian spectrum is in \([0,24]\), e.g. by
\(2|x_ix_j|\le x_i^2+x_j^2\). Thus
\[
 3v^{\mathsf T}L_g^+v\ge \|v\|^2/8.
\]
At diameter squared below \(3/4\), every difference has integral norm
squared below six. No coordinate difference can have magnitude two:
compensating that coordinate requires squared mass at least two on the
others. Translate the cluster coordinatewise by its minima. Every
coordinate is now in \(\{0,1\}\), and every point has the same Hamming
weight. By complement symmetry it is enough to consider weights one,
two, and three (weights zero/six are singletons), with respectively
6, 15, and 20 vertices. In weight three, complementary pairs have
distance squared at least \(3/4\), so the maximum possible size is ten.

There is also a useful small-diameter simplification. At diameter
squared at most \(1/2\), any non-root integer difference would need
norm squared exactly four and an eigenvector at 24. Equality in the
spectral bound requires opposite coordinates across every positive
edge. The graph is connected, so this eigenvector has equal magnitudes
in all six entries. A nonzero integer vector of this form has norm
squared at least six, contradiction. Thus all differences are roots
\(e_i-e_j\) of \(A_5\). Translating one point to zero, two such roots
differ by another root only if they share their positive or their
negative coordinate. In a pairwise family there must be one common
positive coordinate or one common negative coordinate. Up to sign
and translation these are subsets of the six weight-one labels.

[The independent integer checker](verify_clusters.js) evaluates all
55,176 graphs. For each graph and each possible cluster size \(s\), it
computes the greatest integer squared-diameter numerator strictly
beating the incumbent. One then simply checks for a clique of size
\(s\) in each of the relevant constant-weight subgraphs. The radius
argument is exactly \(12-\beta\) with no numerical diameter relaxation;
all density comparisons are integer powers. A floating estimate is
used only to initialize a search for the integer threshold, whose two
bounding comparisons are exact.
