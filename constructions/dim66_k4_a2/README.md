# Dimension 66: 15 translates from a 4-simplex and a triangle

This explicit nonlattice periodic packing has squared minimum **229/30** and
center density

    (15/sqrt(96000)) (229/120)^33.

Its occupied fraction is approximately **2.593018122992473e-13**. It improves
on the immediately previous [16-translate construction](../dim66/README.md)
by **4.876426%**, and on the earlier 9-translate result by **6.659994%**.
The underlying change is a four-simplex times a triangle: despite having
one fewer translate than the previous product, it recovers enough radius
to gain more density.

## Data and exact replay

Use the archived [Nebe–Sloane catalogue Gram matrix](../../research/novel_routes_20261113/higher_parent/neb72.html)
(Git blob 53c8920d3d656821bce869062849001c08137c0e, SHA-256
92daeaf0d93672e82ca5e5bc4479e5775daf59029aa4e983fb464e36ef002044).

    node verify.js ../../research/novel_routes_20261113/higher_parent/neb72.html

The standalone checker verifies all 105 normal-pair distances, the exact
ambient equality witness, primitivity, ambient determinant and source
digest, nonlattice scope, and the exact strict gains. The published
ambient squared minimum 8 is invoked from Gabriele Nebe, *An even
unimodular 72-dimensional lattice of minimum 8*, J. Reine Angew. Math.
673 (2012), 237–247 ([preprint](https://arxiv.org/abs/1008.2862)).

## Period and translates

In the catalogue basis G, take the following six rows of ambient
coefficients A, with zero-based positions:

    (71:+1)
    (59:+1, 66:-1)
    (52:+1, 54:-1)
    (51:+1, 65:+1)
    ( 1:+1,  7:-1, 13:+1)
    ( 5:-1,  9:+1, 13:-1)

The Gram M = A G A^T is the orthogonal direct sum of 10I4−2J4 and
[[8,−4],[−4,8]]. Their Gram determinants are 2000 and 48; det M is
96000. Columns (71,59,52,51,1,5) form the diagonal minor
diag(1,1,1,1,1,−1), so the normal lattice is primitive and the label map
x ↦ AGx maps integer ambient coefficient columns onto Z^6.

Let H=ker(AG) within Z^72. With P the orthogonal projection onto the
66-dimensional perpendicular space in the G-metric, we have
P(x)=x−A^T M^(-1) AGx. The projected normal image of the unimodular
ambient lattice is the dual of the normal lattice, with covolume
1/sqrt(det M); therefore the period H has covolume sqrt(det M).

Let T_r={0,e1,…,e_r} in Z^r and use the 15 labels in T_4×T_2. For
each y, put y_1,…,y_5,−y_6 into the 6 pivot positions above (zero
elsewhere) to obtain b_y. Set x_y=G^(-1)b_y, an integral lift of label y.
The complete center set is the union of H+P(x_y) over those 15 labels.
A period basis can be generated in the same way: for every nonpivot
position j, set b_j=1 and use the six pivot positions to solve Ab=0;
then G^(-1)b gives the corresponding basis vector of H.

In the normal metric M^(-1), every distinct pair in T_4 has squared
distance 1/5, and every distinct pair in T_2 has squared distance 1/6.
Distances add across blocks. The squared diameter is 11/30, with
105 pairs split as 15 at 1/6, 30 at 1/5, and 60 at 11/30.
Every lifted difference between different labels is nonzero in Gamma72,
so its projection loses at most 11/30 from squared norm at least 8.
Consequently distinct translates are genuinely distinct and the squared
minimum is at least 229/30; translations within H have squared minimum
at least 8. For equality, the ambient coefficient vector e24+e25+e40
has norm squared 8 and label (0,1,−1,0,1,−1), which is a difference of
valid labels with squared distance 11/30. Its horizontal projection,
up to a period translation, is one of the actual packing differences.

The map y ↦ P(x_y)+H is a homomorphism with kernel M Z^6:
an integral ambient vector in the real normal span is in A^T Z^6 by
the unit minor. The verifier checks that the class of 2e1 is not one
of the 15 supplied classes, even though the class of e1 is present.
Thus the centers are not a lattice.

## Comparisons

Multiplying the center density by pi^33/33! gives the occupied
fraction above. The square of the ratio to the 16-site predecessor is

    225·90000·229^66 / (256·96000·228^66) > 1.

For the older 9-site construction [verified here](../../research/novel_routes_20261113/higher_parent/gamma72_verified.json),
the corresponding square is

    225·12288·229^66 / (81·96000·225^66) > 1.

The older separately attributed Gamma72/Leech rank-six section derived
from Nebe's Section 4.2 has occupied fraction 2.272528652057e-13;
the improvement over that comparator is approximately 14.102769%.
This result is a constructive lower bound, not an optimality or
lattice-record claim.
