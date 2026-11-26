# A 39-dimensional colour budget and a finite 64-colour target

**No new packing record is claimed.** The comparison is the lattice of [Dorofeev–Sun–Wang, Table 1](https://arxiv.org/html/2607.20359v4#S1.T1), with centre density $3^{16}/(2^{20}\sqrt5)=18.35925948\ldots$. We use the [five-E7 carrier L35](../e7_hermitian35_20261123/) of squared minimum 8, centre density 4, and the [published 38D extension](../e7_lattice38_20261123/). All distances below are Euclidean.

## The obstruction does not depend on the vertical Gram matrix

Let V be any rank-four Euclidean lattice and let f:V→R³⁵/L35 be a homomorphism of finite image size M. Over every site v∈V place the full horizontal translate L35+f(v). The pure vertical lattice K=ker f has index M in V. If the resulting 39D lattice has squared minimum at least 8, then so does K. The exact four-dimensional Hermite constant gives covol(K)≥32 (equality for twice D4). Since covol(K)=M·covol(V), the 39D centre density is at most
$$
  4(\sqrt2)^4/\operatorname{covol}(V)
  =16/\operatorname{covol}(V)\le M/2 .
$$
In particular M≤32 cannot beat even 16, regardless of a deformation of V. The old horizontal colours A,B,C generate a group of order 16. An integral E7⁵ horizontal phase Q is killed by 2 modulo L35, so adjoining it to that group gives at most 32 colours. This argument is for a homomorphic assignment of full L35 translates; it is not an upper bound on arbitrary periodic packings.

## What a 64-colour packing would have to prove

Put D4={z∈ℤ⁴: sum zᵢ is even}, t=√(13/10) and K=2tD4. Enlarge the old group to a subgroup G≅ℤ/4×ℤ/4×ℤ/2×ℤ/2 of the horizontal torus, for example by a class Q of order four whose double is outside the old group. Seek an injective group homomorphism φ:G→D4/4D4. Write
$$
 H(g)=\min_{x\in g+L35}\|x\|^2,\qquad
 \mu(c)=\min_{u\in c+4D4}\|u\|^2 .
$$
For a map φ scaled vertically by t/2, use the union of the 64 cosets (g+L35)×((t/2)φ(g)+K). It is a lattice. Its exact squared minimum is
$$
 \min\{8,\,\lambda_1(K)^2,\,\min_{g\ne0}[H(g)+(13/40)\mu(\phi(g))]\},
$$
since horizontal and vertical minimisations within a coset are independent. Here λ₁(K)²=8t²>8. The 64 translates give covolume covol(L35)·covol(K)/64, with covol(K)=32t⁴. Thus **if and only if** all 63 nonzero-colour costs are at least 8, the prescribed packing has minimum 8 and centre density
$$
 \frac{16\cdot64}{32t^4}=\frac{3200}{169}=18.9349\ldots,
$$
about 3.14% above the DSW benchmark. **The embedding with these 63 costs has not been found.**

## Necessary order-two test

The 15 nonzero elements of G[2] must map bijectively to the nonzero elements of 2D4/4D4. Writing such an element as 2y, its cost μ is 4 times the squared coset minimum of y mod 2D4. The 24 roots of D4 occur in 12 opposite pairs in different cosets, giving twelve costs of 8. The remaining three cosets are represented by 2e₁, (1,1,1,1), and (−1,1,1,1), which all have squared minimum 4: none is zero or contains a root, and their displayed representatives have squared norm 4. Those three and zero form a plane under addition, giving precisely three costs of 16. Consequently, every order-two colour with H(g)<27/5 must belong to a common plane in G[2], and H(g)≥14/5 is mandatory for all order-two colours. After that inexpensive filter, the full 63-colour embedding remains to be tested.

## The surviving doubled syndrome

The integral part of the old colour group is $I=\{0,911,4096,5007\}$. The old colour 5007 has horizontal squared depth 4 and therefore takes one of the three exceptional order-two slots. An integral colour is cheaply certified to have depth at most 4 if it is zero, one nonzero nonradical E7 phase, or two root phases in distinct blocks. The checker below generates 7,939 such syndromes. For each coset $d+I$, where $d=2Q$, the binary span of 5007 and its cheap syndromes must fit in the exceptional plane. This leaves 18 cosets, each with exactly two cheap colours whose sum is 5007. In particular the preliminary example $d=117$ is excluded: all four integral colours in its coset are cheap.

Each remaining plane is already full. The colour $B+d$ is another order-two colour outside it. A finite integer recurrence finds a representative of $2(B+d)$ with squared norm 16 for 17 of the 18 remaining representatives, proving horizontal squared depth at most 4 and excluding them. The only doubled-syndrome coset not excluded is $1496+I$, i.e. 1496, 1623, 5592 or 5719. This gives a necessary search filter, not a proof of a lower bound on the depth of the survivor or of a 64-colour embedding. It also holds for other common scales $t$ in this model that could beat DSW: $t^2\ge2$ is already too sparse, so an ordinary order-two slot gives a total strictly below 8 when horizontal depth is at most 4.

Run `node verify_double_filter.js` for the integer recurrence and witnesses. It uses the published E7 and 38D certificates. Its finite-box costs are **upper** witnesses; no exhaustion of the whole horizontal lattice is required for the exclusions.

## Scope of the direct 38→39 searches

For a uniform translate q of the actual L38 at layer spacing h the centre density is 16/h. In a general periodic four-layer parallelogram 0,q,q+r,r with alternating gaps a,b, the centre density is 32/(a+b). If T(s)=sqrt(max(0,8−dist(s,L38)²)), the smallest possible a+b is max{T(q)+T(r),T(q+r),T(q−r),√2}. In either approach at least one adjacent horizontal torus displacement needs squared clearance >7.2404967349 to beat DSW.

An exact min-plus distance oracle calibrated on the 16 cosets of the 38D certificate screened 16,000 half-coordinate and 20,000 quarter-coordinate horizontal patterns of the five E7 blocks, each with 8,192 integral syndrome options and 128 half-grid positions in the 3D component. Best observed adjacent squared clearance was 7 on this grid. A separate box-cover search on 12,000 half-coordinate patterns also failed to reach the threshold when the 3D component ranged over the full torus. Its cell upper bound at midpoint m and half-width r, using a particular periodic target image p, is horizontal squared cost + ||m−p||²+2r||m−p||₁+3r²; minimise over the 16 colours and periodic images. The reported box certification concerns only those chosen horizontal masks. There is no all-phase or all-lattice claim here.

Run the archived [search replay](replay.md) for the scoped experiments.
