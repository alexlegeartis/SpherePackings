# Multiple phases above each E7 site: the same 19D bound

The [single-phase theorem](README.md) allows at most one horizontal K12 coset at each projected site. The **same all-radius ceiling** `0.08449535563928...` holds if each site instead carries any subset of the dual palette. Keep the same hypotheses otherwise: normal sites `sqrt(s) E7` with scalar metric, dual phase labels, and an untwisted `2 sqrt(s) E7` period. This still does not rule out label-shifted or larger periods, other metrics, or off-palette phases.

Let `mu` be squared separation, `M` the total number of site/phase pairs in the period and `F(mu,s)` be the density factor defined in the main proof. Thus `delta=(M/128)F(mu,s)` and every nonempty packing has `s≥mu/8`. For `mu>8/3` the largest nonzero dual-coset depth is too shallow for two phases at one site, so the main proof applies. We treat `mu≤8/3` here.

## The nontrivial window: 2 < mu ≤ 8/3 and s < mu/6

Every phase occurs at most once in the projected period, since all normal residue classes have representatives of squared norm at most 6. At a fixed site, all nonzero phase differences must have `q=1`. Across root pairs, `mu−2s>(2/3)mu>4/3`, so `q=2` is forbidden. Use the ordered horizontal relations `(zero, nonzero q=0, q=1, q=2)` and normal relations `(zero,R,N,T)` of least norms `(0,2,4,6)`. The only allowed nonidentity product relations, possibly including more than the actual geometry allows, are

```
(1,1),(1,2),(1,3), (2,0),(2,1),(2,2),(2,3), (3,2),(3,3).
```

Take the unnormalized character eigenmatrices of the two relation schemes:

```
P = [[1,224,252,252],  Q = [[1,63,63, 1],
     [1,  8,-18,  9],       [1,-1,-1, 1],
     [1,  8,  9,-18],       [1, 9,-9,-1],
     [1,-19,  9,  9]]       [1,-7, 7,-1]]
```

Rows are frequency classes with a common tuple of relation sums. If a group element is in shell `h`, the average of its character over horizontal frequency row `i` is `P[i,h]/P[0,h]`; analogously for normal row `j` and shell `g`. The [checker](verify_multiphase.py) verifies both this shell-wise uniformity and the eigenmatrices directly, without reliance on a numerical SDP or a claim about automorphism orbits.

For each relation let `X[h,g]` count ordered pairs of packing layers with that difference, divided by `M`. These entries are nonnegative, `X[0,0]=1`, and `sum X=M`. Each product frequency row average of the squared absolute character sum of the layer set, divided by `M`, is nonnegative. Hence, with `c_ij(h,g)=−P[i,h]Q[j,g]/(P[0,h]Q[0,g])`,

```
sum_(allowed nonidentity h,g) c_ij(h,g) X[h,g] ≤ 1.
```

Use the rows and nonnegative rational weights

```
(i,j):    (0,2)    (1,0)    (1,2)    (3,0)    (3,2)
weight:       7    126/5    882/5    112/5    784/5
```

Their combined coefficient is at least one at each of the nine allowed positions. This gives `M≤1+1939/5=1944/5`, hence `M≤388`. There is no additivity or random-independence assumption on any packet of labels. Consequently in this window

```
delta ≤ (388/128)(8/27)(2/3)^6 = 1552/19683.
```

## Other low-radius windows

* If `2<mu≤8/3` and `s≥mu/6`, each site is a clique for horizontal `q=1`, of size at most seven. Thus `delta≤7 (sqrt(3)/16)(2/3)^6`.
* If `4/3<mu≤2` and `s≥mu/6`, each site's phase set is independent in the `q=2` graph. It has degree 252 and least eigenvalue −18 (read the last column of `P`); the eigenvalue argument of the main proof with zero edges gives cap `floor(729·18/(252+18))=48`. Thus `delta≤48 (sqrt(3)/16)(1/2)^6=3 sqrt(3)/64`. If instead `s<mu/6`, every phase is used at most once, so `M≤729` and `delta≤(729/128)(8/27)(1/2)^6=27/1024`.
* If `0<mu≤4/3` and `s<mu/2`, a given phase's support in the normal quotient is root-independent. It uses at most 12 sites by the main theorem's normal least eigenvalue, so `M≤729·12` and `delta≤(729·12/128)(8/27)(1/3)^6=1/36`. If `s≥mu/2`, use `M≤729·128` and `F(mu,mu/2)=F(mu,mu/8)/128`, giving `delta≤1/432`.

The largest low-radius bound is `3 sqrt(3)/64≈0.08119`, strictly below the stated single-phase all-radius ceiling. The exact comparison and dual certificate replay with `python verify_multiphase.py`.
