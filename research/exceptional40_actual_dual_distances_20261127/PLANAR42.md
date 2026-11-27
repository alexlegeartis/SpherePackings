# Exact 42D planar-stack optimum over the exceptional 40D carrier

**No new density record.** Let `L` be the rank-40 orthogonal section of the [actual-distance theorem](README.md), with squared minimum 6, covolume 81 and finite phase group `A=L*/L`. Take any union of full translates of `L` in `R^40 × R^2`, with horizontal phases in `A` and arbitrary real planar sites. Several phases may use the same site if the packing constraints allow it. There is no additivity, periodicity or fixed-radius hypothesis. The center density is at most

```
3^(35/2) / 2^21.
```

The bound is attained by the `sqrt(3) E6` section of [Dorofeev–Sun–Wang, Theorem 1, §§2.1, 3.1 and Appendix A.1](https://arxiv.org/html/2607.20359v4). Their nested primitive `sqrt(3) E6 ⊂ sqrt(3) E8 ⊂ P48p` sections ensure that the 42D lattice is a union of full translates of `L`; its horizontal projections are in `L*` by integrality. The result also holds for any exceptional 40D carrier satisfying the preceding actual-distance theorem, even if a particular carrier is not presented together with that 42D equality construction. This is not an upper bound for arbitrary 42D packings or for real horizontal phases outside `L*/L`.

## 1. Distance caps and density normalization

Write `mu` for the squared packing separation, with `0 < mu ≤ 6`, and `I` for upper planar site intensity, counting each horizontal phase at a site as a separate layer. A pair of layers of phase difference `u` and planar separation `d` has minimum squared interlayer distance `H(u)+d²`, where `H(u)` is the actual horizontal squared coset depth. Here `H(0)=0`. Any open disk of radius `sqrt(mu−T)/2`, for `0 ≤ T < mu`, contains only pair differences with `H(u)>T`, and hence distinct phases. The previous theorem gives these caps on its population:

| `T` | cap `M` | reason |
|---:|---:|---|
| `0` | `6561` | all phases |
| `10/3` | `2187` | avoid the two nonzero elements in one bad cyclic direction |
| `4` | `81` | avoid nonzero isotropic differences |
| `14/3` | `9` | every difference has normal finite quadratic value 1 |

For the final row, let `q` denote the `F3`-valued normal quadratic form `3q_K` and `B` its polar form. Translate a phase to 0. The other `j` phases have polar Gram `I_j+J_j`, diagonal 2 and off-diagonal 1, because their individual and difference quadratic values are all 1. For `j=9` this Gram is invertible in dimension eight, impossible. We also need its exact equality structure: **any 9-clique of `q=1` differences has phase sum zero.** For `j=8` the Gram has rank seven, with kernel the all-ones line. The eight phase vectors cannot span the full nondegenerate eight-dimensional space, since then their Gram would be invertible; they span dimension seven by the rank lower bound. Thus their coefficient kernel equals the Gram kernel and their sum is zero. Restoring the translated phase does not affect the sum, as `9=0` in `F3`. In particular two distinct 9-cliques cannot share eight phases.

Averaging the disk count yields `I ≤ 4M/(pi(mu−T))`. This and the later weighted disk inequality are valid for nonperiodic sites: integrate over the planar square enlarged by the disk radius, so every layer in the original square contributes a whole test disk area; then take the upper limit as the square grows. The resulting bound is uniform in the center of the square. On each horizontal translate the lattice count in a horizontal cube is its volume divided by 81 plus a uniform boundary error. Hence the full upper point intensity is at most `I/81`, first in product cubes and then in balls or other standard averaging regions by cube subdivision. Since the radius is `sqrt(mu)/2`, the ratio of center density to the value in the theorem is at most

```
G(mu) I/sqrt(3),     G(mu)=(mu/6)^21.
```

The four disk bounds give `D_T(mu)=4M G(mu)/(pi sqrt(3)(mu−T))`. A separate pair-spacing bound for `mu>16/3` follows by writing `a=mu−16/3`: all planar sites are distinct and separated by at least `sqrt(a)`. Thue's planar circle-packing density theorem yields `I ≤ 2/(sqrt(3) a)` and ratio

```
P(mu)=2 G(mu)/(3 a).
```

The intensity version of Thue applies to arbitrary sites: periodically repeat any finite square of sites with an extra border of width `sqrt(a)` to obtain a periodic packing, then take the upper limit over arbitrarily large squares.

## 2. Crossing-annulus lemma

The pair and disk bounds alone overlap slightly above one. The following geometric argument closes exactly that gap. Assume `mu>16/3`, rescale planar distances by `1/sqrt(a)` so their minimum is 1, and put

```
R = sqrt((mu−14/3)/(4a)),       e=37/1000,       S=R+e.
```

Assume `31/20 ≤ R ≤ 33/20`. Inside an open test disk of radius `R`, all phase pairs have `q=1` because their actual horizontal depth exceeds `14/3`. The same is true of any pair, inside or outside, whose normalized separation is strictly less than `2R`. All sites within radius `S` have **distinct phase labels**: equal labels need unscaled separation at least `sqrt(mu)`, whereas `4 S² a ≤ 4(1687/1000)²(2/25)<16/3<mu`. The bound `a<2/25` follows by solving for `mu` from `R≥31/20`.

Call an inside–outside pair **long** if its normalized distance is at least `2R`. Here inside means radius `<R`, and outside means the annulus from `R` to `S`. These cross pairs form a matching. To prove this, fix a site of radius `r` and let the opposite sites have radii at most `s`. To be long an opposite site must belong to the circular cap cut off by projection at least

```
h(r,s) = (4R²−r²−s²)/(2r)
```

onto the direction opposite the fixed site; there are no long pairs if `r=0`. For `s≤S`, `h(r,s)` decreases with `r`. Thus the two relevant cap heights are bounded below by `h(S,R)` and `h(R,S)`, which are positive. The cap is contained in a disk of squared radius `s²−h²` centered at the foot of this cut. Put `U=3R²−S²`. The larger of the two squared cap diameter bounds is

```
4S²−U²/R² = 16 R e + 4 e² − 4 e³/R − e^4/R²
            < 16(33/20)e + 4e² < 1.
```

The other squared bound is `R²/S²` times this one. Since distinct sites are 1-separated, each vertex of the cross graph has at most one long partner.

Each endpoint of a long cross pair lies in the ring `[R−e,R+e]`, by the triangle inequality. This ring contains at most ten sites total. Otherwise two among eleven, ordered around its center, have angular separation at most `2 pi/11`, and their squared distance is at most

```
4e² + 4 S² sin²(pi/11) < 4e² + 16(1687/1000)²/49 < 1,
```

using `pi<22/7`. All annulus sites are in this same ring.

Let `f` count sites inside radius `R`, and `g` count sites in the annulus. Then

```
5 f + g ≤ 45.
```

Indeed `f≤9`. If `f=9`, an outside site would be short to all inside sites but possibly one. It would create either a 10-clique or a second 9-clique sharing eight phases with the existing one. Thus `g=0`. If `f=8`, classify outside sites as matched (`m`) or unmatched (`u`); then `u≤1` because two unmatched sites would each join the eight inside sites in a 9-clique. The matching endpoints and any unmatched site all lie in the ring, so `2m+u≤10`; with `u≤1` this gives `g=m+u≤5`. If `f≤7`, use the ring cap `g≤10`. Equalities on the radial boundaries cause no problem (they also have zero area in the averaging).

Undoing the normal rescaling and averaging this weighted count gives the sharpened `T=14/3` disk ratio

```
D_(14/3)(mu) / C(R),      C(R)=1+((1+e/R)^2−1)/5.
```

Namely the weighted effective area is `pi a [R²+(S²−R²)/5]`, instead of `pi a R²`. This improvement uses the equality structure of the finite phase clique; the annular count is not just a bound for an uncolored multiple disk packing.

## 3. Exact all-radius finish

Set `mu(R)=16/3+(2/3)/(4R²−1)`, and use the following rational cuts:

```
muL = mu(33/20),       mu1 = 10817/2000,       mu0 = 540883/100000.
```

On successive intervals ending at `15/4`, `41/10`, `19/4` and `muL`, use the four basic disk caps from the table in order. For fixed `T>0`, the logarithmic derivative of `G(mu)/(mu−T)` is `21/mu−1/(mu−T)`. It can change sign only from negative to positive, so the maximum over a closed subinterval is at an endpoint. For `T=0` the bound increases from zero. Each needed positive endpoint is strictly below one by rational comparisons using lower bounds for `pi` and `sqrt(3)`.

The entire interval from `muL` to `mu0` satisfies `31/20≤R≤33/20`; the basic `T=14/3` disk ratio is increasing there. On `[muL,mu1]` use `C(R)≥C(33/20)>D_(14/3)(mu1)`. On `[mu1,mu0]`, the decreasing `R` satisfies `R(mu1)<1571/1000`, and use `C(R)≥C(1571/1000)>D_(14/3)(mu0)`. The cuts are chosen for the **shifted** crossing: the spacing bound itself was not strengthened.

Finally on `[mu0,6]` use `P(mu)`, whose same logarithmic-derivative argument gives an endpoint maximum. Exactly `P(mu0)<1` and `P(6)=1`. This proves the bound, strictly if `mu<6`. The cited DSW cross-section has determinant `3^7` and squared minimum 6, giving `(6/4)^21 / 3^(7/2)=3^(35/2)/2^21` as claimed.

Run `python verify_planar42.py` for a standard-library, rational-only replay of the Gram kernels, crossing-annulus constants, exact cuts and endpoint signs. There is no arrangement enumeration or numerical optimum assumption.
