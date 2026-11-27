# Sphere157: bare cross-sections cannot improve by keeping a mixed quaternionic tail

**No new density record.** The statement controls the explicit 51D–53D quaternionic parents at the inherited packing radius, along **all rational section directions**, including those oblique to their vertical coordinates. It does not cap the actual packing radius of a section missing every norm-six vector.

## Setup and statement

Let P be any even unimodular 48-dimensional lattice with minimum squared norm 6. Let V be the root lattice D_k (minimum squared norm 2 and Gram determinant 4), 3≤k≤5. For any additive homomorphism from V to the ambient horizontal space modulo P, the resulting **graph extension** is a lattice L ⊂ R^48⊥R^k with full kernel P and projection V:
```
 L = {(p+f(v), v): p∈P, v∈V}.
```
Here a lift f modulo P can be chosen linearly on a basis of V. Assume L has minimum squared norm 6, as for the exact [project 51D–53D constructions](https://github.com/alexlegeartis/SpherePackings/tree/main/research/p48n_hurwitz_20261122): their three leading vertical Gram blocks are D₃, D₄, D₅ (change signs on a tree of simple roots), their horizontal determinant is one, and their Gram Schur determinant is 4.

Throughout, det denotes the **Gram determinant**, i.e. squared covolume. A full rank-d lattice from a bare section has center density (3/2)^(d/2)/sqrt(det) at radius sqrt(6)/2. Define H_μ(j) as the infimum of Gram determinants of rank-j lattices of minimum squared norm at least μ, with H_μ(0)=1.

For a rank-d rational section B of L, take its primitive dual normal K=L*∩span(B)^⊥. Let r be the rank of the purely vertical kernel of its horizontal projection; let s=k−r. Thus s measures the number of normal-space vertical directions **not** removed by purely vertical dual constraints; for an oblique section this is not assumed to equal the dimension of B's intersection with the vertical coordinate space. The horizontal projection of K has rank j=48−d+s. The determinant factorization yields
```
             det B ≥ H₆(48−d+s) · H₂(s).                 (1)
```
It does not depend on the phase labels f. If s=0 the section is an ordinary bare cross-section of P, not an outcome of mixed gluing. No section at inherited radius that uses s>0 beats the named density benchmarks in dimensions 39–48. In 40–48 even s=0 never improves those benchmarks. Dimension 39 with s=0 remains the base problem of a rank-nine normal to P with determinant below 21870.

For dimensions 49–53, the sharp inherited-radius optimal determinants in this parent family, whenever its dimension permits that section, are respectively 2, 3, 4, 4, 4. In particular no oblique normal with nonzero horizontal projection can improve the project’s entries there.

## Proof of the determinant factor

The standard primitive dual-section identity is det B = det L · det K. Here is a basis proof. Since B is a saturated rational section, extend a basis of B to a lattice basis of L. Orthogonally projecting the complementary basis vectors gives a basis of the projected lattices in span(B)^⊥. The projected lattice's dual in that plane is exactly L*∩span(B)^⊥: pairing a vector of that plane with any L vector is pairing with its projection. The Euclidean covolume of L is the product of covolumes of B and the projected lattice. Squaring these identities proves the claim. This also proves the analogous identity below for V*.

The projection of L* to the horizontal space lies in P*=P because L contains (P,0); its purely vertical sublattice is exactly (0,V*) because the vertical projection of L is all of V. In particular the horizontal projection K_p of K is a discrete rank-j sublattice of P of squared minimum at least 6 if j>0. The vertical kernel K_v is a saturated rank-r section of V*: it is V* intersected with an actual rational vector subspace. Its orthogonal complementary section in V has rank k−r=s and minimum at least 2 if s>0. The primitive dual-section identity gives
```
 det K_v = det(V*) · det(V∩span(K_v)^⊥) ≥ H₂(s)/4,
```
where the last orthogonal complement is taken **inside the vertical space**. This equality also holds for r=0 or s=0 using determinant one in rank zero.

Choose a basis for K_v and lifts in K of a basis for K_p. Projecting those lifts off span(K_v) shows that their Gram is G_p+E, where G_p is the projected horizontal Gram, and E is the positive-semidefinite Gram of the remaining vertical components. (A coordinate Gram blocks argument gives the same Schur factor.) Hence det K ≥ det K_v det K_p; this is an inequality even for arbitrarily tilted graphs of lifts. Since det L=det P det V=4, the factor 1/4 cancels and proves (1).

When s=0, K_v spans the entire vertical space, so its section is inside P×{0}; the only remaining orthogonality constraints are to K_p. When the horizontal projection has rank j=0, K is purely vertical. Then B is the graph extension over the remaining section of V and includes the whole P. Its true minimum is 6 because it inherits the parent minimum and also contains P.

## Numerical thresholds stated as exact inequalities

We only need standard low-dimensional Hermite inputs and one published low-dimensional packing bound:
- The exact squared determinant minima at squared minimum two in ranks 0 through 8 are (1,2,3,4,4,4,3,2,1), realized by the standard root lattices. These are the classical exact Hermite constants (see the [Nebe–Sloane catalogue of densest lattices](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/density.html)). Therefore H₆(j)=3^j H₂(j) for j≤8.
- The primary [Cohn–Elkies 2003 exact-arithmetic LP bounds, Table 3](https://annals.math.princeton.edu/wp-content/uploads/annals-v157-n2-p09.pdf#page=23) give **center** density <3/50 in rank 10 (tabulated about 0.05804) and <2/25 in each of ranks 9 and 11–14. These deliberately loose rational values are used directly as center densities, with no unit-ball-volume approximation. For min squared 6 this gives H₆(10)≥(3/2)^10(50/3)^2=4100625/256>16018 and H₆(j)≥(3/2)^j(25/2)^2 otherwise in 9–14.

For target density δ_d define the determinant threshold T_d=(3/2)^d/δ_d², so improvement requires **det B<T_d** at the inherited radius. The following exact threshold formulas (or strict integer upper bounds) suffice. For 39–43 the source is [Dorofeev–Sun–Wang, v4, Table 1](https://arxiv.org/html/2607.20359v4#S1); for 44,45,47 use [Chen–Hu–Li–Wang–Wu](https://arxiv.org/html/2505.02394v1#S3); for 46 use [Conway–Sloane, The Antipode Construction](https://neilsloane.com/doc/anti.pdf); 48 is P itself.

| d | Threshold T_d or a proved upper value | smallest lower value from (1) over s=1,…,5 (downward integer bound) |
|---:|---:|---:|
| 39 | 21870 | 32036 |
| 40 | 6561 | 12013 |
| 41 | 4374 | 13122 |
| 42 | 2187 | 8748 |
| 43 | 2^180·3^85/137^43 <728 | 4374 |
| 44 | (3/2)^44·605·(110/157)^44/25 <216 | 1944 |
| 45 | 2^138·(3/23)^45 <55 | 648 |
| 46 | 3(27/26)^46 <18 | 216 |
| 47 | 2^189·(3/47)^47 <6 | 54 |
| 48 | 1 | 12 |

The s=0 bound at d40–42 equals the comparator, and at d43–47 its values 972,324,108,27,6 respectively exceed the thresholds. For d48 it is exactly 1. For d39 the proof does not settle j=9 and s=0; (1) is taken only as its exact reduction to the already-understood bare P section problem. This distinction matters.

For d=49,…,53, j=s−(d−48)≥0. The cases j=0 have the exact determinant lower values H₂(d−48)=2,3,4,4,4; they are attained by vertical sections of types A₁,A₂,D₃,D₄,D₅ as their ranks permit. For A₁ use a primitive root, A₂ the first three coordinates with sum zero, and for D_j use the coordinate subspace; these are primitive. The base P is intact, so the min6 radius is attained. For j>0, H₆(j)≥6 and H₂(s)≥3 in the relevant ranges, giving determinant ≥18, strictly worse. The live [49D archive](https://github.com/alexlegeartis/SpherePackings/blob/main/research/covering_and_deep_cosets_20261115/deeper_results.json), [50D exact proof](https://github.com/alexlegeartis/SpherePackings/blob/main/research/adaptive_layers_20261115/P48_PROOF.md), and [51D–53D proof](https://github.com/alexlegeartis/SpherePackings/tree/main/research/p48n_hurwitz_20261122) have those comparator determinants.

**Scope.** These are bare rational cross-sections only; higher actual section minimum, projected packs, additional antipode cosets, nonlattice modifications, other vertical lattices, dilations coupling the two orthogonal spaces, and the bare rank-nine normal problem in 39D are not excluded.

## Exact arithmetic replay

Run `node verify.cjs`. The standalone standard-library script checks all 50 positive-tail comparisons, the zero-tail comparisons, and the tight pure-vertical targets, using only exact BigInt rational arithmetic. Hosted replay of the same source passed. The determinant factor itself is a basis identity, not an enumerative claim.