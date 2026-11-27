# Uniform 27D gap: this planar family is below 87.6% of Vardy

This independently strengthens the conclusion in [the main proof](README.md)
for precisely the same carrier. There is still no new packing record.
Use all of its notation, in particular \(A_0=\mu-28\), \(t=\mu/32\),
and the eight retained fractions \(p_f\) of the heavier bank. Set
\(c=219/250\), and suppose the actual ratio \(R\) to the incumbent is at
least \(c\). The period-area estimate and main density formula give
\[
 W\ge\frac{cV}{4t^{27/2}}\ge4c\sqrt3\,t^{-25/2}>6,
 \qquad t^{25}\ge 3c^2/4>(39/40)^{25}.
\]
On this range the logarithmic derivative calculation in the main proof
again gives \(A_0\ge4t^{27/2}\): indeed
\(378-(25/2)(32\cdot39/40)<0\).
In particular, if \(V\ge A_0 M\), then \(W\ge cM\).

## Caps detect the Hamming-three constraints

For any set of \(k\) labels with no bad edges, rescale the transverse
distances by \(A_0^{-1/2}\). The main Voronoi bound gave area at least one
at each site, with up to four full caps of area \(g=\pi/8-1/4\).
Every positive-area cap neighbor is of the other colour and is at
distance strictly below \(\sqrt2\). Two different translates of a
single label cannot both cut caps of one site: their mutual squared
distance would be below \(8\) in normalized units, whereas the period
minimum is at least \(\mu/A_0\ge8\), as \(\mu\le32\).

Thus if a site has at most \(4-e\) possible opposite labels that can
cut caps, its normalized cell area is at least \(1+eg\), simply by
bounding each cap by \(g\). In counting possible opposite labels one
should exclude any complementary three-bit label \(f+\mathbf1\) with
\(p_f+p_{f+\mathbf1}>1\): its matched horizontal difference has weight
five and forces squared normal distance at least \(\mu-20>2A_0\).
Call such a retained pair of complementary labels *safe* for the cap
count. We will use \(g>57/400\), which follows from
\(\pi>157/50\) (e.g. the classical lower estimate \(223/71\)).

## The graph and weighted deletion

Write \(\ell_f=1-p_f\), total loss \(L=8-W<2\). Two disjoint bad edges
are impossible since each needs total endpoint loss at least one.
A pairwise intersecting family of edges is a star or the three edges
of a triangle, except for the empty graph.

* **Triangle.** Delete two triangle labels and keep the other six; there
  is no bad edge remaining. If both endpoints of a surviving
  complementary pair are outside the triangle, their losses sum to
  less than one because the entire triangle has loss at least \(3/2\).
  If one endpoint is the remaining triangle label, the two deleted
  triangle labels again consume at least one and the pair's loss is
  less than one. At least two of the original four complementary pairs
  survive, so at least four cells lose one possible cap on top of the
  missing caps from the colour counts. The colour counts are \(3+3\)
  or \(4+2\); their missing-cap totals before counting complements are
  respectively six or eight. Hence \(V\ge A_0(6+10g)\). This gives
  \(W\ge c(6+10g)>13/2\), whereas the triangle requires \(L\ge3/2\).
* **Nonempty star.** Retain the seven labels other than the star
  center. Their counts \(4+3\) force at least four missing caps, so
  \(W\ge c(7+4g)>13/2\) and \(L<3/2\). For a single bad edge choose
  as center the endpoint of loss at least \(1/2\), which leaves loss
  below one among the seven. If the star has at least two edges,
  the center's loss is at least \(2-L>1/2\) from the two edge
  inequalities, with the same conclusion. The remaining seven labels
  include three whole complementary pairs and all are safe. The six
  additional missing caps give \(V\ge A_0(7+10g)\), so
  \(W\ge c(7+10g)>7\), inconsistent with any bad edge.
* **No edges.** The ordinary area bound yields \(W\ge8c>7\).
  All four complementary pairs are safe and every cell has at least
  one missing cap. Hence \(V\ge8A_0(1+g)\); the density formula now
  gives \(R\le1/(1+g)<c\), a contradiction.

For a compact exact replay of the numerical signs used above:

~~~python
from fractions import Fraction as Q
c = Q(219, 250)
g0 = Q(57, 400)
assert 48*c*c > 36
assert 3*c*c/4 > Q(39, 40)**25
assert 378-Q(25, 2)*32*Q(39, 40) < 0
assert c*(6+10*g0) > Q(13, 2)
assert c*(7+4*g0) > Q(13, 2)
assert c*(7+10*g0) > 7
assert 8*c > 7
assert 1/(1+g0) < c
~~~
