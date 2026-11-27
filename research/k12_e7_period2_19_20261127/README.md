# A spectral bound for 19D Coxeter–Todd / E7 colored stacks

**No new density record.** We tested a potential 19-dimensional construction: stack full translates of the Coxeter–Todd lattice `K12` over a scaled E7 lattice, choosing which projected sites to use and assigning each a dual-coset phase. At the most direct scale and radius, using all sites would give center density `sqrt(3)/16`, about 21.7% above the current benchmark. Deleting one site per eight still leaves headroom. But for the compact period described here, two finite spectra jointly rule it out, including all scalar normal scales and all radii.

## Family and comparator

Fix `K12` with squared minimum 4 and covolume 27. In the orthogonal normal space use sites `sqrt(s) E7`, with `s>0` and E7 at squared minimum 2. At each site one may omit the layer, or include all of `K12+a` for a phase `a` in the finite group `A=K12*/K12`. Require that the occupied sites and their phase labels repeat under `2 sqrt(s) E7`, without a horizontal label shift. The period then has 128 projected residue classes. There is at most one phase per site. Placements need not be additive or use all sites.

The center density in this family is at most

```
(74/128) (8/27) (8/9)^6 = 0.08449535563928...
```

This is a ceiling, not a claim that the value is attained. The external [Chen–Hu–Li–Wang–Wu 19D antipode construction](https://arxiv.org/html/2505.02394v1#S3) has center density `13^(19/2)/(3^9 5^(21/2)) = 0.08895810212922...`, corresponding to the [pinned project row](https://github.com/alexlegeartis/SpherePackings/blob/3a8e374c74497e5efaaeea7cb46d17ff1588b530/LOWER_BOUNDS_1_96.md). Even the older [laminated lattice `Lambda19`](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/LAMBDA19.html) has center density `1/(8 sqrt(2)) = 0.08838834764832...`. At `s=2/3` and squared separation 4 specifically, the bound is stronger: at most 82 of the 128 sites can be used, whereas at least 106 are required to beat the Chen et al. benchmark. We do not exclude larger or label-shifted periods, different normal metrics, real phases outside the dual palette, multiple layers at one site, or other carriers.

## Two finite graphs

The Coxeter–Todd horizontal phase quadratic form is `q(a)=3||a||²/2 mod 3`. The lattice is even 3-modular (in particular `sqrt(3) K12*` is isometric to `K12`), and its squared covering radius is `8/3`; see [Conway–Sloane, Corollary 7 and the dual description after (44)](https://neilsloane.com/doc/Me107.pdf). Therefore `A` is nondegenerate of dimension six over `F3`. Its Gauss sum divided by `sqrt(|A|)` is `-1` by the even-lattice signature identity, since the signature is 12. The nondegenerate form with polar Gram `I_6+J_6` has the same sign (its exact value counts are `(225,252,252)`), hence is isometric to it by the two-type classification over `F3`. In this model

```
q(v) = 2 (sum_i v_i² + (sum_i v_i)²) mod 3.
```

Because a nonzero dual vector has squared norm at least `4/3`, and its norm is at most `8/3` in a closest representative modulo `K12`, the phase depths are exact: zero for the zero class, 2 for nonzero isotropic classes, `8/3` for `q=1`, and `4/3` for `q=2`.

The horizontal graph joins two palette labels if `q(a−b)=1`. Its degree is 252, and its character eigenvalues, with multiplicity, are `252:1, 9:476, −18:252`. For the character with frequency `w`, put `z_w = #{v : q(v)=1, w·v=0}`. The eigenvalue is `(3 z_w−252)/2`, since the connection set is invariant under negation. The accompanying checker evaluates all 729 dot-product counts exactly. This graph has clique number at most seven: translating a vertex to zero, seven other vertices of a hypothetical 8-clique would have polar Gram `I_7+J_7`, invertible over `F3`, impossible in dimension six.

For a coset model of the projected graph, write `E7` as the sum-zero vectors in `R^8` all of whose coordinates are integral or all half-integral. It is the index-two extension of `A7` by `c=(1,1,1,1,−1,−1,−1,−1)/2`, and has covolume `sqrt(2)`. A basis is `c,a_2,...,a_7`, where `a_j=e_(j−1)−e_j` in zero-indexed coordinates. (The first standard `A7` coefficient of `c` is `1/2`.) Thus `E7/2E7` is a seven-bit vector space, with coordinates given by these integer basis coefficients mod 2.

The 126 norm-two vectors give 63 distinct nonzero labels `R`, and the 756 norm-four vectors give 63 different labels `N`. The last nonzero label `t` is the label of `(3,3,−1,−1,−1,−1,−1,−1)/2`, of norm six. Exact coordinate enumeration below gives `R+t=N`. Since `E7` is even, these are the actual least norms in the four cases. We call `{x,x+t}` a twin pair. Between any two distinct twin pairs, the root edges are exactly one perfect matching: the two possible differences are `r` and `r+t`, one in `R` and one in `N`.

The normal graph joins projected residue classes if their difference is in `R`. It is 63-regular and has spectrum `63:1, 9:28, −1:63, −7:36`. These are the 128 exact binary character sums `sum_(r in R) (−1)^(w·r)`.

## The transfer inequality

Let `M` be the number of occupied projected classes, `e` the number of induced normal root edges, and `m_a` the multiplicity of the horizontal phase `a`. From the least normal eigenvalue,

```
2e >= (70/128) M² − 7M.                       (1)
```

Indeed, decompose the occupancy indicator into its constant projection and a vector perpendicular to constants, and apply the adjacency matrix's least eigenvalue on the latter. In particular if root edges cannot occur, then `M≤12`.

Suppose now that every occupied root edge requires a horizontal `q=1` difference, and norm-four normal pairs cannot share a phase. Then the horizontal graph must have at least `e` edges after counting phase multiplicities, where each phase pair contributes the product of its two multiplicities. The largest nonconstant horizontal eigenvalue bounds that count from above:

```
2e <= M²/3 + 9 sum_a m_a².                     (2)
```

To see the constant `1/3`, use `(252−9)/729`. No label occurs more than twice, and any repetition is on a twin pair. Two distinct repeated labels lie on different twin pairs, which are joined by root edges, so the repeated labels form a horizontal `q=1` clique. Thus there are at most seven such labels and `sum_a m_a²≤M+14`. Combining (1) and (2) gives `41 M² /192 ≤16 M+126`. Every integer `83≤M≤128` violates it, so `M≤82`. If norm-six pairs also cannot share a phase, all labels are distinct and `41 M²/192≤16M`, giving `M≤74`.

## All squared separations and scales

Write `mu` for the required squared separation, `0<mu≤4`. The same-phase normal period has squared minimum `8s`, so a nonempty packing requires `s≥mu/8`. A pair with unscaled squared normal length `j` needs horizontal depth at least `mu−js`. The center density is

```
delta = (M/128) F(mu,s),  F(mu,s) = (mu/4)^(19/2) / (27 sqrt(2) s^(7/2)).
```

In particular `F(mu,mu/8)=(8/27)(mu/4)^6` and `F(mu,mu/6)=(sqrt(3)/16)(mu/4)^6`.

* If `mu−2s>8/3`, no root edges are possible. With `M≤12` and `s≥mu/8`, `delta≤1/36`.
* If `2<mu−2s≤8/3`, occupied root edges force `q=1`, and `4s<mu` because `mu≤4`. If `6s≥mu`, use `M≤82` and `s≥mu/6`; this gives `(82/128)sqrt(3)/16`. If `6s<mu`, use `M≤74` and `s≥max(mu/8,(mu−8/3)/2)`. These two scale bounds cross at `mu=32/9`. The density bound increases up to there, and then decreases through `mu=4`: on the second scale branch the logarithmic derivative of its square is `19/mu−7/(mu−8/3)`, negative up to `38/9`. The largest value is `(74/128)(8/27)(8/9)^6`.
* If `mu−2s≤2`, use `M≤128` and `s≥max(mu/8,(mu−2)/2)`. The first branch increases to `mu=8/3`. On the second branch the square has logarithmic derivative `19/mu−7/(mu−2)`, with just a minimum at `19/6`. Thus the global bound in this case is the larger of `512/19683` and `1/(27 sqrt(2))`.

The largest of all these ceilings is the stated rational number, strictly less than the Chen et al. benchmark by exact squared comparison. The [standard-library checker](verify.py) replays the coset enumeration, both graph spectra, the occupancy inequalities and all endpoint comparisons.
