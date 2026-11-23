# A quantitative 39D obstruction for the five-E7 carrier

This is a family bound, not a bound on arbitrary packings. Use a horizontal copy of the [35D lattice L](../e7_hermitian35_20261123/) with exact squared minimum 8. Its horizontal phases are integral shifts in E7^5 plus any choice of the half-radical in each of the five blocks. Place one full translate of L at every site of any four-dimensional Euclidean lattice V. The phase assignment may be nonlinear and may have a period larger than V. At squared minimum 8 its **39D center density is less than 16.36**. This is below the prior center density 3^16/(2^20 sqrt(5)), about 18.35926.

The proof also permits moving every vertical site by Euclidean distance at most 0.004 and moving each horizontal phase by distance at most 0.005 in the horizontal torus. The center-density cap then becomes **18.32**, still below the prior bound. Larger motions, other fractional phases, partial vertical site sets, and minimum below 8 are not excluded.

## Exact horizontal depth table

A phase difference is described by five fractional radical bits and a 13-bit integral syndrome. For fractional Hamming weights j=0,1,2,3,4,5, its squared horizontal torus depth is at most respectively 6, 5.5, 5, 6.5, 6, 7.5. The precise integer min-plus computation checks all 262,144 differences. `verify_phase.py` imports the hash-pinned [E7 form and alphabet](../e7_hermitian35_20261123/certificate.json), checks an inverse-Gram coordinate bound, and compiles a complete local squared-norm-14 enumeration. It then verifies the quarter-unit binary syndrome convolution for every fractional pattern. In particular, the weight-three maximum 6.5 occurs on exactly 72 of 8192 integral syndromes for every triple.

## Short vertical directions

Write M(j) for the caps just listed. We work with the distances in the underlying vertical lattice V, even in the perturbed statement. The following inequality holds for the unperturbed vertical squared separation q of any two sites whose assigned fractional patterns differ in j entries:
```
q >= 8 - M(j) - eta,
```
where eta=0 for no motions and eta=11/100 for the stated motions. The horizontal distance to a torus period is 1-Lipschitz. For motions of size epsilon, its squared upper bound therefore increases by at most 11*epsilon+4*epsilon^2 since sqrt(M(j))<=11/4. A vertical motion of size zeta at both endpoints can reduce the squared-separation requirement on V by at most 12*zeta: subtract 2*zeta in length, square, and use sqrt(8)<3. The two losses sum to at most eta. A bounded site displacement does not change the intensity of its periodic arrangement, and separation remains positive.

The shortest nonzero squared length a of V is at least 1/2-eta. If a>=3/2-eta, the exact four-dimensional lattice Hermite bound makes the vertical covolume at least a^2/2, so the center density is at most 32/a^2. Otherwise every edge in the shortest direction v must change its fractional pattern by all five radical bits, the only difference with M(j)>6.5. In particular v is primitive. Project V perpendicular to v, producing a three-dimensional lattice T. The vertical covolume is sqrt(a) times the covolume of T; the center density of the 39D packing is 16 divided by this vertical covolume.

Each projected site has a well-defined fractional color modulo complement, because its preimages along v alternate between two complementary patterns. If there is only one color modulo complement, globally only two patterns are used. Each fixed pattern's subset of V has squared separation at least 2-eta, even before displacement. Using the published repository's rigorous general **4D occupied-density bound <0.63**, its vertical point intensity is at most 16 c/(2-eta)^2 where c=2*0.63/pi^2. The two colors combined give a 39D center density at most 512 c/(2-eta)^2. This case also proves the no-motion two-pattern exclusion for arbitrary periodic vertical sites, not only a full lattice. Separately, if all vertical squared distances are at least 1.35, the same 4D bound without color separation rules out any improvement.

## The three-dimensional parity-cut certificate

For a nonzero projected difference of squared length s, its shortest lift without a parity restriction has longitudinal squared length at most a/4. Thus s>=3a/4. To use the fractional colors, choose the parity of lifts with the smaller Hamming weight j<=2. Its shortest longitudinal length is t in [0,sqrt(a)]; the other parity's shortest longitudinal length is sqrt(a)-t. Therefore
```
s >= max(8-M(j)-eta-t^2, 8-M(5-j)-eta-(sqrt(a)-t)^2).
```
Minimizing over t gives 2-a-eta for j=0, 3-a-eta for j=2, and 9/4-a/4-1/(16a)-eta for j=1 when a>=1/2 (2.5-a-eta otherwise). These inequalities give a floor l=max(3a/4,2-a-eta) for all projected vectors, and a higher floor h=min(3-a-eta, j1_bound) for differences between distinct projected colors.

By Selling's theorem T has an obtuse superbase of four vectors. Write its six nonnegative conorms as p_ij. For each of the seven nonzero parity classes, the minimum squared norm is the associated cut sum: sums of superbase vectors attain it, because the norm of a general lift with coefficients n_i is the sum of p_ij(n_i-n_j)^2, with crossing differences odd. A parity-cut vector shorter than h preserves the projected color under *every* translation. If their labels spanned the three parity bits, take three independent cuts with determinant odd. Three-dimensional 0/1 determinants have absolute value at most two; hence those three cuts form an integer basis and every site has the same color. In the nonconstant-color case there must therefore be a nonzero binary functional vanishing on all short cuts. Its four nonzero-value classes have norm at least h. This reasoning does not assume any additivity of the phase assignment.

Enumerate the seven choices of functional. For each, vertices of the six-conorm polyhedron use six independent equalities selected from six nonnegativity and seven cut lower constraints. Its recession cone is nonnegative. The squared covolume is the spanning-tree polynomial of the conorm graph, so it is monotone along recession directions. Concavity of log determinant on positive definite Gram matrices implies a lower bound attained at one of the finite constraint vertices. A rational partition of the a range into intervals of length 1/1000 supplies uniform lower rational bounds on l and h. `verify_selling.py` checks exact lower bounds on the squared covolume of V of 0.988455983303 without motion and 0.8176613445 with eta=11/100, in the respective small-a, nonconstant-color cases. Squaring the comparisons gives the stated center-density caps without floating-point acceptance.

## Replay and inputs

Run `python -B verify_phase.py` (C++17) and `python -B verify_selling.py` (NumPy). The arithmetic is exact. The separate [4D certificate](../../verifications/dim4/refinement_energy20261115/) supplies the general 4D upper bound. The fixed baseline is the [39D entry in the pinned general catalog](https://github.com/alexlegeartis/SpherePackings/blob/ab3ce28ea92d3b265502de6e4c2c02bc3f6fe7ce/LOWER_BOUNDS_1_96.md), from the [2026 cross-section and antipode paper](https://arxiv.org/abs/2607.20359). Selling's theorem and the conorm model are in [Conway–Sloane, Voronoi reduction of three-dimensional lattices](https://neilsloane.com/doc/Me169.pdf).
