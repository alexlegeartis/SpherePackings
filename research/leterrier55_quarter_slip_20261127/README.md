# 55D: a noncollinear shear gate and exact blockers on four root orbits

**No new density lower bound.** We tried to beat the verified [three-layer 55D packing](../uneven_layers55_20261122/) by moving the long-gap phase away from the old order-three line. A single missing nearest-contact pairing would have yielded a **76/75-fold** improvement, from center density 92,197.2234 to 93,426.5197. The large initial contact sample suggested two possible orbits. Exact polynomial replay gives a blocker for each, and in fact rules out the whole one-slip family on four disjoint root orbits at every radius. It does not show that other slip directions, periods or carriers cannot improve.

The 54D parent `L` has minimum 10 and determinant `3^25` by [Leterrier, Corollary 1.4](https://doi.org/10.1007/s40993-022-00321-0). The two nonzero order-three cosets `±Q+L` have actual squared depth `D=28/3`; use the [validated source-to-Gram binding](../cross_coset55_20261116/PROOF.md), not the faulty 55D lamination. The current period height is `H0=(19/5)h`, with `h²=2/3`, and it has three layers. The [thesis's printed 116,783.1496](../tilted55_20261116/PROOF.md#external-source-correction) is an unverified comparison: its displayed construction has an explicit norm-six vector. Neither the hypothetical primal gain nor the dual gain below exceeds that printed number.

## The candidate and the precise contact bottleneck

For a norm-10 vector `v in L`, use horizontal increments `Q,Q,Q+t v`, vertical gaps `h,h,H−2h`, and period `(t v,H)` modulo `L` (`3Q∈L`). At `t=1/4` set `H=(15/4)h`. The period vector has squared norm 10. For a closest displacement `w in Q+L`, the pairing `w·v` is an integer because `Q∈L*`, and `||w+v||²≥D` implies `w·v≥−5`. If `−5` did not occur, every closest `w` would satisfy

```
||w+v/4||² ≥ D−2+10/16 = 191/24.
```

The long vertical gap has norm `49/24`, so the third one-gap comparison would meet 10. Every nonclosest `w` has squared norm at least `D+2` because `L` is even and `Q∈L*`, and interpolation between `w` and `w+v` gives an even larger bound `D+3/2−15/8`. For two gaps, the unmodified pair gives `D+4h²=12`; either modified pair is bounded by `D−15/8+(11h/4)²=25/2` since both endpoints of the horizontal interpolation belong to `−Q+L`. Three gaps give the period vector and at least 10, since `||v/4||≤sqrt(10)/2` forces zero closest to the fractional period shift. Four or more gaps give at least `16h²>10` using heights alone. Thus the one contact-avoidance condition would give the stated exact `H0/H=76/75` gain.

Conversely one closest `w` with `w·v=−5` makes the exact horizontal distance profile

```
dist²(Q+t v,L) = D−10t(1−t)     (0≤t≤1).
```

The lower bound is convex norm interpolation over each pair `w,w+v`, both in the same coset; that contact realizes equality. It also supplies the corresponding contact for the negative slip direction by replacing `w` with `w+v`. In particular the proposed quarter shear has a literal squared distance `19/2`, not 10, across its long gap.

## Exact root-orbit contacts from the curve

The parent is the narrow Mordell–Weil lattice for `y²=x³+x+t^28` over `F_(3^6)(t)`. Use [the published point arrays and basis data](../cm58_20261116/frontier/) and the [existing complete Q-contact enumeration](../uneven_layers55_20261122/LOCAL_OPTIMUM.md). The coefficient field is represented modulo `a^6+2a^4+a²+2a+2`, with coefficient digits in base three. The four exact narrow norm-10 polynomial candidates are minus point #28 and points #30, #4, #0 (zero-indexed); their archived root rows start at 0, 2187, 4374 and 6561, respectively. All have `deg x=10, deg y=15`. In this model a polynomial point with `deg x≤9` and `deg y=14`, leading `y` coefficient `+1`, has height `28/3` and is in `Q+L`. Indeed the two nonidentity components at infinity are distinguished by the signs of the leading `y` coefficient (in the chart `x_s=s^10 x(1/s), y_s=s^15 y(1/s)`, they have initial `y_s=±s`). This is also checked by the earlier Gram binding.

Here is a cheaper exact test for the needed pairing. Let `v` be a narrow norm-10 polynomial point and `w` a closest polynomial point of `Q+L`. Put `d=x_v−x_w`, `e=y_v−y_w`. Then `deg d=10`, and the rational x-coordinate of the sum has numerator `N=e²−(x_v+x_w)d²` and denominator `d²`. The leading degree-30 terms cancel by the curve equation, so `deg N≤29`. The common factor of `N,d²` is the common factor of `e²,d²`, of degree `2 deg gcd(d,e)`. The height formula therefore gives `height(v+w)=D+20−2 deg gcd(d,e)` and hence

```
w·v = 5−deg gcd(d,e).
```

In particular the exact narrow-side test is simple polynomial divisibility: `d` divides `e` if and only if the pairing is `−5`. Both points are on the curve, both coset assignments are pinned above, and at divisibility `w+v` is again precisely a norm-`D` point in `Q+L`.

We use all 2,187 isometries `t→t+a` with the three constant lifts of each correction `R(t)`, where `R³+R=(t+a)^28−t^28`. These lifts form a group: composition adds translation parameters and any difference of corrections is one of the three constant solutions of `R³+R=0`. They preserve the narrow lattice and `Q+L` by the infinity signs. The [verifier](verify.mjs) independently validates the literal blocker for each representative (and its inverse transport), and checks that the four root orbits are disjoint of size 2,187, thus giving **8,748 distinct blocked directions**. For the two promising directions the 2,187-point initial subshell missed every blocker. The optional full replay uses the exhaustive 1,143 normalized closest-point seeds and pairs them with all 2,187 transforms of each candidate; it finds exactly **341, 258, 341, 260** minus-five contacts in the four cases. The earlier full Q enumeration justifies the normalized coefficient equations and its completeness; this verifier directly checks the field arithmetic, seed count and literal obstruction. No claim is made about all 15,309,000 roots of `L`.

## No radius can rescue this one-slip geometry on these orbits

Keep the increments `Q,Q,Q+t v` but allow arbitrary nonnegative gaps at any squared separation `0<mu≤10`. Modulo `L` and changing `v` to `−v` as needed, take `0≤t≤1/2`. Write `H` for the period height. The 2D period sublattice generated by `(v,0)` and `(t v,H)` has minimum at least `mu`, so planar lattice packing gives `H≥sqrt(3) mu/(2sqrt(10))`. Against the live construction, a gain requires `H<H0 (mu/10)^(55/2)`. For `mu≤99/10` the squared ratio from the planar bound is at most `(H0²/(15/2))(99/100)^53<1`, so such a gain is impossible.

For `99/10≤mu≤10` set `u=mu/10` and `g=sqrt(mu−D)`. The period vector requires `10 t²≥mu−H²`. Under the alleged gain this implies `t≥tB=sqrt(14/375)`, because `mu−H0² u^55` is decreasing on this interval, with value `28/75` at `mu=10`. In fact its derivative is at most `1−(55/10)·9·(46/100)<0` by Bernoulli. The two short gaps then each require length at least `g`; the last gap requires length at least `sqrt(g²+10t(1−t))` by the exact contact profile. Since `t(1−t)` increases on `[0,1/2]`,

```
H ≥ T(mu) = 2g + sqrt(g²+10 tB(1−tB)).
```

At `mu=10`, the excess squared term is greater than `112/75`, because `tB>14/75`; thus `T(10)>H0`. For smaller `mu` in this interval the derivative of `T` is at most `3/(2g)<2` since `g²≥17/30`. The derivative of the required height `H0 u^(55/2)` is greater than `3·(55/20)·(73/100)>6`, again by Bernoulli. Therefore `T(mu)>H0 u^(55/2)` throughout. This excludes every radius and gap choice on those four root orbits, not just the particular test at `t=1/4`.

## A distinct positive gate left open

If some shortest dual vector `R` of norm `10/3` satisfies `dist²(Q+R,L)≥26/3`, use increments `Q,Q,Q+R/2`, gaps `h,h, sqrt(55/6)−2h`, and period `(R/2,sqrt(55/6))`. Horizontal norm interpolation gives at least `49/6` in the long gap's coset; together they give `20−(4/3)sqrt(55)>10`. The modified two-gap horizontal distance is at least `11/2`, since `2Q+R` is nonzero and the [verified dual minimum](../leterrier54_dual_normal_20261126/) is `10/3`; vertical addition gives `46/3−(2/3)sqrt(55)>10`. The other gaps are checked as above. This would improve the **verified** comparator by factor `sqrt((722/75)/(55/6))≈1.0248`. The root-orbit certificates do not address this dual condition, and no qualifying dual vector is asserted.

Replay quickly with `node verify.mjs`; for the full contact counts use `node verify.mjs --full`. Neither path needs FLINT, a closest-vector solver, or floating-point distance tests.
