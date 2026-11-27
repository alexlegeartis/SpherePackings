# Exact 26D lattice optimum for the Bacher-22 carrier with any normal metric

**No new density result.** Let `S` be the particular rank-22 section of Bacher26, of covolume 8 and squared minimum 4, described in the [22+4 decomposition](../2026-11-26-27d-bacher-glue/README.md). Among **all lattices** in 26D containing `S` horizontally and with horizontal coordinates in `S*`, the greatest center density is **1/√3**, attained by Bacher26. This allows every four-dimensional normal period metric, skew between the period directions, every subgroup of horizontal classes, and any sphere radius. A global similarity can first match the given scaled copy of `S` to this normalization. The result concerns this one carrier, not all 26D lattices or nonlinear periodic packings in 26D.

The same horizontal palette in **28D** cannot beat 2/3 even with an arbitrary periodic normal stock at each label. In **27D** there is a separate obstruction to a projected short `A3` root packet, but an improvement over Bacher within the entire 27D lattice family is not proved impossible.

## Inputs, exact geometry, and separate benchmarks

The source 26D Gram is [Bacher (1996)](https://doi.org/10.1007/s002220050180) in [Martinet's literal basis collection](https://jamartin.perso.math.cnrs.fr/Lattices/bacher.gp). Use the exact 64-class phase table in the [22+4 decomposition and verifier](../2026-11-26-27d-bacher-glue/). The accompanying [local-rigidity proof](../bacher_local_rigidity_20261126/) also verifies a norm-4 vector in `S`. Write `H=S*/S=Z4²×Z2²`, with class `a+4b+16c+32d` in the earlier generator order. For nonzero `h` the least horizontal squared distance of this class is `4−W_h/4`, where `W_h∈{2,3,4,6,7,8}` is the table value (we use `W` for the **integer** table in this note). Thus any lattice of squared minimum `mu=4−lambda` in this carrier must lift each nonzero class `h` to squared normal distances at least `W_h/4−lambda`. For zero label, a nonzero normal displacement must separately have squared length at least `mu`.

The live and attributed center-density comparisons are: **26D** Bacher26 `1/√3` for both lattice and unrestricted tables; **27D** Bacher lattice `1/√3`, Vardy periodic unrestricted `1/√2`; **28D** Bacher lattice `2/3`, Vardy periodic unrestricted `1`. See the [Nebe–Sloane catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/density.html), [Cohn's table](https://cohn.mit.edu/sphere-packing/#data), [Vardy (1999)](https://doi.org/10.1090/S0002-9947-99-02169-8), and the [live project unrestricted table](../../LOWER_BOUNDS_1_96.md). At a fixed dimension all comparisons of center density and occupied fraction have the same ratio. The Bacher-family 27D tie with `sqrt(2) A5` is [already certified](../bacher_oblique_and_a5_20261126/), and gives `1/√3`, not the unrestricted 27D comparator.

## Fiber bound (including low radius and 28D periodic)

Normalize `S` to squared minimum 4; since it is a horizontal period, `mu≤4`. Fix a horizontal label. All points with this label have horizontal coordinates in one coset of `S`, so by applying periods every pair can be compared at the **same horizontal point**. Thus the corresponding normal stock must itself pack at minimum squared `mu`. Its density in the full packing is the normal point intensity divided by `covol(S)=8`. More explicitly, a full period group may be taken to contain `S`; passing to finite index fixes each of the finitely many horizontal labels, and projecting that subgroup to the normal space gives a full period group for the stocks there. Its density factor is justified by the block triangular basis of the period with a horizontal basis first. For a lattice the stocks for all labels in use are translates of the same normal lattice: they are fibers of the homomorphism from the packing lattice to `H`. Consequently, with `m` labels and normal dimension `r`, the center density is at most

```
(m/8) (mu/4)^11 delta_r
```

where `delta_r` is the corresponding unrestricted or lattice normal bound. The exponent 11 is the horizontal packing-radius factor; the normal radius has already been included in `delta_r`. For 26D lattices `delta_4=1/8` (equivalently `gamma_4=sqrt(2)`). With `m≤64` and `mu≤15/4`, this gives `(15/16)^11<1/√3`. If `m<64`, then the lattice image in `H` is a subgroup, so `m≤32` and the same inequality gives `delta≤1/2` for **every** radius.

In 28D use the rigorous 6D **center-density** upper bound less than `0.081` from [Cohn–Elkies, Table 3](https://annals.math.princeton.edu/wp-content/uploads/annals-v157-n2-p09.pdf), which displays `0.08084`. It gives the full arbitrary-periodic center-density ceiling `8·0.081=0.648<2/3`. In 27D the **lattice** normal bound `delta5=1/(8 sqrt(2))` only gives `1/√2`, so it rules out beating Vardy but **does not** rule out beating Bacher. The unrestricted normal bound in 5D is different; we do not substitute a lattice bound for it.

## 26D high radius: the entire additive palette, now as a projected lattice

It remains to assume `15/4<mu≤4` and full image `H`. Project the packing lattice onto the normal space to get `Q`. Two distinct horizontal labels cannot occupy the same projection since their horizontal squared distance would be at most `7/2<mu` (choose a translate of `S`). Hence each projection has a unique label, giving a well-defined surjective homomorphism `ell:Q→H`. Any nonzero projected vector has squared length at least `1/2−lambda>0`, so `Q` is discrete and full rank four. The horizontal kernel consists exactly of `S`: otherwise compare the extra residue to zero. Lifting a basis of `Q` proves `covol(M)=8 sqrt(det Q)`. The induced map `Q/2Q→H/2H` is a surjection between binary rank-four spaces, hence an isomorphism. This step is specific to the four-dimensional normal here, not the 27D case.

Let `tau=3/4−lambda`. Any nonzero vector of `Q` of norm (squared length) `<tau` must have `W=2`. These labels occupy exactly three nonzero parities in `H/2H`: `2,6,10`, when the four low-order bits of `a,b,c,d` encode the parity. They are linearly independent. If two such projected vectors `v,w` have the same parity and are not antipodal, both `(v+w)/2` and `(v−w)/2` are nonzero vectors of `Q`, and the average of their norms is `(||v||²+||w||²)/4<tau/2<1/2−lambda` since `lambda<1/4`. This contradicts the minimum of `Q`. Thus there are `r≤3` short antipodal pairs, represented by vectors `v_i` independent modulo 2 and therefore over the reals. If the two short parities differ, the palette shows their sums and differences have `W≥4`, hence norms at least `1−lambda`. It follows that, with `a_i=||v_i||²`, their Gram has off-diagonal entries of absolute value at most `(a_i+a_j−1+lambda)/2`.

Let `A0=Σ Z v_i` in their real span. The squared covering radius of `A0` within that span is at most `Σ a_i/4`: round the coordinates of any point independently, rounding up with probability equal to its fractional part, and the mean squared rounding error has zero cross terms and diagonal contribution at most this bound. This also proves primitivity in `Q`. An extra vector of `Q` in the span but not in `A0` would produce, after rounding, a new nonzero vector of norm at most `Σa_i/4<tau`, not among the listed pairs. The orthogonal projection of `Q` modulo this primitive subspace is a lattice of rank `4−r` (for `r=0` use `Q` itself); lifting and rounding any nonzero projection gives squared normal minimum at least

```
f = tau − (sum_i a_i)/4 > 0,
```

because that lifted nonzero vector is outside the span of all short pairs. Projection factors the Gram determinant of `Q` as the product of the determinant on `A0` and the Gram determinant of the complementary projected lattice.

## A small exact determinant inequality

Use `a_i=1/2−lambda+x_i`, `0≤x_i≤1/4`, `x_i+x_j≥lambda`. Let `D(x)` have these diagonal entries and off-diagonal entries `−(x_i+x_j−lambda)/2`. This matrix is strictly diagonally dominant with positive diagonal (also at the closed endpoints). In rank two the Gram determinant is at least `det D`, by the off-diagonal bound. In rank three, expand the determinant: the three negative square terms are minimized at maximum absolute values, and the product term is bounded below by minus twice the product of those maxima. Hence `det D` is again a lower bound. For rank zero use determinant 1. The Hermite constants give

```
det Q >= det D(x) · f^(4−r) / gamma_(4−r)^(4−r),
1/gamma_(4−r)^(4−r) = (1/4,1/2,3/4,1) for r=0,1,2,3.
```

For fixed `lambda`, the logarithm of the right side is concave on the indicated polytope: `log det` of a positive affine matrix and `log` of a positive affine function are concave. Its minimum is thus at a vertex. Up to permutation, the vertices for each rank are:

```
r=0: ()
r=1: (0), (u)
r=2: (0,lambda), (0,u), (u,u)
r=3: (0,lambda,lambda), (0,lambda,u), (0,u,u), (u,u,u), (lambda/2,lambda/2,lambda/2)
u=1/4
```

A vertex with zero coordinate has each remaining coordinate at one of its resulting interval endpoints. With a coordinate at `u`, the remaining lower-dimensional constraints give the listed cases. Without coordinate bounds active in rank three, all three pair constraints must bind, giving the last vertex. This list remains valid when constraints coincide at an endpoint.

At each of the eleven vertices compare the rational polynomial above with `(3/64)(1−lambda/4)^26`. Set `z=4 lambda∈[0,1]`. For the difference polynomial in the power basis `d_j z^j`, its Bernstein coefficients of degree 26 are `b_k=Σ_(j≤k) d_j binom(k,j)/binom(26,j)`. Every such coefficient is nonnegative by the exact standard-library checker; all with `k≥1` are strictly positive. Since the Bernstein basis is nonnegative on the unit interval, this proves `det Q≥(3/64)(mu/4)^26`. The full density is `(mu/4)^13/(8 sqrt(det Q))≤1/√3`. The known Bacher26 lattice, comprising all 64 glued classes in the original period, attains equality with `mu=4`. This proves the entire asserted 26D lattice scope.

## 27D: a short `A3` packet obstruction and the actual scope of the attempted gain

At `mu=4`, suppose the projected lattice contains three vectors `b1,b2,b3` with Gram `t(I_3+J_3)`, where `t<7/16`. These vectors form with zero a tetrahedron of squared side `2t<1`, so the six differences must all have nonzero labels with `W≤3`. **These four-cliques do exist** in the palette: there are exactly 144 unordered triples with zero as the fourth vertex. However each has a forbidden extra short combination. If `a,b,c` are the labels of `b1,b2,b3`, at least one of `a+b−c`, `a+c−b`, `b+c−a` is zero or has `W≥7`, as the exact finite replay checks. In the `A3` Gram the corresponding vectors are nonzero and have squared norm `4t<7/4`, less than the required value `4` for label zero or `W/4≥7/4` otherwise. This excludes the whole short root packet, not all short tetrahedra. A projected rank-five homothetic `A5` or `D5` normal is then impossible at the lattice-improvement gate `det Q<3/64`: its scale `t` would have `det Q=6t^5` or `4t^5` respectively, forcing `t<7/16`; the three differences from a fourth common coordinate give the specified `A3` Gram. Arbitrary rank-five normals are not excluded.

The same finite graph has no five-clique. Consequently no five projected points, even in a nonlinear periodic packing over the same palette at squared separation 4, can have all their squared pair distances below 1. The additional `A3` packet constraint above uses additivity and is stronger for lattice constructions.

We also tried projection with three orthogonal half-norm coordinates, labels `4,20,36`, and a two-dimensional oblique tail of Gram `[[3/8,1/8],[1/8,1]]`, with horizontal lift `(1/2,1/2,1/2)` and label 1 for the first tail basis vector, zero horizontal lift and label 2 for the second. Its determinant would give a **1/24 reduction compared to the Bacher-family determinant**, only about **2.15% density gain** over `1/√3` if valid, still below Vardy. It is **not valid**: subtract twice the first tail vector from the second, a zero-label vector of squared length `7/2<4`. Neither this one lift failure nor the short-packet obstruction rules out general 27D lattice or periodic improvements.

Run `python verify.py` for all palette identities, eleven Bernstein certificates, comparison arithmetic and the explicit failed tail witness. No numerical inequality, native enumeration or new density record is used.
