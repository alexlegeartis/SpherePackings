# Dimension 66: 16 translates from two orthogonal Gamma72 normals

We construct an explicit nonlattice periodic sphere packing in dimension 66
with minimum squared distance **38/5** and center density

    (16/300) (19/10)^33.

Its occupied fraction is approximately **2.472450874704875e-13**,
a strict **1.700637%** increase over the previously verified 66-dimensional
9-translate construction (2.431106576447802e-13). This is a constructive
lower bound, not a claim of optimality or a new lattice bound.

## Certificate and sources

The ambient lattice is Gabriele Nebe's even unimodular Gamma72 of minimum
squared distance 8, from *An even unimodular 72-dimensional lattice of minimum
8*, J. Reine Angew. Math. 673 (2012), 237–247
([author's preprint](https://arxiv.org/abs/1008.2862)).
Use the complete, row-major 72 by 72 Gram matrix G in
[Nebe–Sloane's catalogue](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/neb72.html).
The archived source is
[neb72.html](../../research/novel_routes_20261113/higher_parent/neb72.html)
(Git blob 53c8920d3d656821bce869062849001c08137c0e, SHA-256
92daeaf0d93672e82ca5e5bc4479e5775daf59029aa4e983fb464e36ef002044).

    node verify.js ../../research/novel_routes_20261113/higher_parent/neb72.html

The replay checks raw input provenance, symmetry, even integrality, positive
definiteness and unimodularity, the literal normal witness and saturation,
all 120 projected-normal distances, a literal distance-equality witness,
nonlattice scope, and the exact strict density comparison. It relies on
Nebe's published minimum-8 theorem for the ambient lattice.

## Construction

Use integer coefficient columns x for the basis of Gamma72, with squared
norm x^T G x. Let A have the following rows, where unlisted coefficients
are zero and positions are zero-based:

    (50:+1, 62:-1)
    ( 3:+1, 62:+1)
    ( 0:+1, 16:+1)
    (25:+1, 43:-1, 45:-1)
    (57:+1, 69:+1, 71:+1)
    ( 1:-1, 21:+1, 27:-1, 45:+1)

Write M = A G A^T. Direct multiplication yields M=K⊕K, where

    K = [[ 8,-3,-3],
         [-3, 8,-2],
         [-3,-2, 8]],       det K = 300.

The columns (50,3,0,25,57,1), in that order, are
diag(1,1,1,1,1,-1), so the six-row normal is primitive and A maps Z^72
onto Z^6. Unimodularity of G implies that the label map x ↦ AGx is
surjective as well. Let E=ker(AG) in R^72 and H=ker(AG) in Z^72.
The normal space is the column span of A^T. Orthogonal projection
(in metric G) onto E is

    P(x) = x - A^T M^(-1) A G x.

For any y in Z^6, take b_y supported on the six specified columns,
putting y_1,…,y_5,-y_6 there in order, and set x_y = G^(-1)b_y. Then x_y
is integral and A G x_y=y. This gives a completely specified period and
translation rule with no numerical coordinate search. An explicit basis
for H can be assembled by taking each nonpivot column j, putting 1 in
position j of a vector b, putting the unique values on the pivot columns
that satisfy Ab=0, and taking G^(-1)b.

The image of Gamma72 under projection to the normal space is the
dual of the six-row normal lattice, via y ↦ A^T M^(-1)y. Its covolume
is 1/sqrt(det M). In orthonormal horizontal and vertical coordinates,
a basis of H followed by any lifts of a normal-image basis is block
triangular. Since Gamma72 has covolume 1, H has covolume
sqrt(det M)=300.

Let T={0,e1,e2,e3} in Z^3 and use the 16 labels in T×T. The packing
centers are the union of H+P(x_y) over those labels. In one block,
the exact normal squared distances between distinct labels in T are
one of the values 11/60 and 1/5, with four at 11/60 and two
at 1/5. In two blocks their squared distances add, so the diameter square
is 2/5. For any two centers in distinct supplied translates, represent
their difference as the projection of an integral ambient vector with
label difference y-y'. It is necessarily nonzero before projection. Its
horizontal squared norm is at least 8−2/5 = 38/5. This also proves
the translates are distinct. Within one translate the minimum is at
least 8. Equality is attained between distinct translates: the ambient
vector with coefficient column e_0+e_12 has squared norm 8, label
(-1,0,0,0,1,-1) in (T×T)−(T×T), and projected normal squared norm 2/5.
Its projection, up to a period translation, is an actual packing
difference. Thus the packing radius is sqrt(38/5)/2.

For clarity, y ↦ P(x_y)+H is a group homomorphism on all of Z^6 with
kernel M Z^6. Indeed, if P(x_y) belongs to H, subtract it from x_y:
the resulting integral ambient vector lies in the real normal span,
whose integral vectors are exactly A^T Z^6 because the displayed minor
is a unit. The converse follows by taking a vector in A^T Z^6.
The replay checks that twice the class of the first unit label is not
in the supplied 16 cosets, so their union is not a lattice.

## Exact comparison

The previous verified 66D construction is the 9-translate Gamma72
certificate gamma72_root_cluster_rank6.npz, listed in
[its verification ledger](../../research/novel_routes_20261113/higher_parent/gamma72_verified.json),
with normal det=12288 and minimum squared at least 15/2. Its reported
center density is (9/sqrt(12288))(15/8)^33. The square of the new-to-old
ratio is the rational number

    (256·12288·76^66)/(81·90000·75^66) > 1.

Multiplying the center density by the radius-one 66-ball volume
pi^33/33! gives the stated occupied fraction. The older, separately
attributed Gamma72/Leech-block rank-six section, [derived from
Nebe's Section 4.2](../../research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py),
has occupied fraction approximately 2.272528652057e-13; the present
construction is about 8.797347% denser than that comparator.
