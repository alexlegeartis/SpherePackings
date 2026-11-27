# Four cusps constrain 54D Mordell–Weil glue beyond integral lattices

No new packing is claimed. Among **all lattices** `L ⊂ M ⊂ L*`, including nonintegral ones, with squared minimum at least 8, a strict improvement over Leterrier's 54D packing must have one of two forms:

| Index `[M:L]` | Squared minimum | Extra requirement |
|---:|---:|---|
| 3 | 10 | even integral; the new coset is isotropic, with exactly 15,707,034 vectors of squared norm 10 |
| 729 | 8 | even integral; the scaled dual has at most 6,480 vectors of squared norm 10 |

Neither existence claim is proved. This includes nonintegral extensions of arbitrarily large index in the stated minimum range, not a screen of sampled candidate vectors. Minimum below 8, over-lattices not contained in `L*`, and different carriers are outside this classification. For **arbitrary nonlinear periodic motifs** in `L*/L` at the original squared separation 10, there are at most **three shifts**, and more than one is possible exactly when the first lattice row is possible. The second lattice row's additional restrictions and a companion per-coset result are below.

## Source and normalization

Take the integral even Leterrier lattice `L` of rank 54, squared minimum 10 and Gram determinant `3^25`; see [Leterrier (2022), Corollary 1.4](https://doi.org/10.1007/s40993-022-00321-0) and his [thesis](https://infoscience.epfl.ch/record/307330/files/EPFL_TH10417.pdf), Theorem 3.2.7 and Proposition 3.3.10. These are the published 54D inputs, not the disputed 55D lamination. The [existing exact Gram and scaled-dual proof](../leterrier54_dual_normal_20261126/) checks that `3G^-1` is integral even, and proves that `J=sqrt(3)L*` has squared minimum 10 and 1,347,192 shortest vectors using both cusp series at level 3. We use these proved inputs. The live and named external 54D comparator are the same Leterrier packing; its center density is `(5/2)^27 / 3^(25/2)` and the [occupied density](../../LOWER_BOUNDS_1_96.md) is about `1.46699394710872e-10` after multiplying by `pi^27/27!`. An index `3^s` lattice of squared minimum `mu` has density ratio `3^s (mu/10)^27` to it, in either convention.

## Deriving the four-cusp identities

Write `q=exp(2 pi i z)`, `Q=q^(1/3)`, `A=L*/L`, with finite quadratic form `q_A(u)=||u||^2/2 mod Z` and polar form `b_A`. Every class is killed by 3 because `3G^-1` is integral, and `q_A` takes values in `(1/3) Z/Z` because `3G^-1` is even. Thus `A` is a nondegenerate quadratic space over `F_3` of dimension 25, with form `3 q_A`. For `u in A`, let `theta_u(z)=sum_{x in L+u} q^(||x||^2/2)`. Its first possible nonconstant `Q` exponent is at least 5, by the scaled-dual minimum. For any real weights `w_u`, write `f=sum w_u theta_u` and let `n_r=sum_{3 q_A(u)=r} w_u`. Initially assume just `f=c+A8 Q^12+O(Q^13)`; for the packing applications this follows from separation at least 8 applied to every positively weighted class. Here `c=w_0`.

The theta-vector is a weight-27 holomorphic modular form for the level-three Weil representation, whose restriction to `Gamma(3)` is trivial; see [Scheithauer, *The Weil representation and some applications*, Theorem 4.7 and the triviality on the principal congruence subgroup used there](https://www.mathematik.tu-darmstadt.de/media/algebra/homepages/scheithauer/publications/Weil_representation.pdf). Here are the needed normalizations directly: the standard `T` multiplies `theta_u` by `exp(2 pi i q_A(u))`; under weight-27 slash by `S`, Poisson gives `(-i)^27 3^(-25/2) sum_v e^(-2 pi i b_A(u,v)) theta_v`. Using the opposite Fourier sign has no effect on any cusp constant or positivity below. In particular at each of the other three cusps represented by `T^k S`, `k=0,1,2`, the `Q^1,...,Q^4` coefficients vanish. This statement does **not** assume `M` is even or integral; the theta components used are those of the even base `L`.

Let `x=theta_(A2)` and `y=theta_(A2+v)` with `v=(1/3,1/3)` in a basis of Gram `[[2,1],[1,2]]`. They are weight-one forms on `Gamma(3)`, `x=1+6Q^3+O(Q^9)`, `y=3Q+3Q^4+O(Q^7)`. For `omega=exp(2 pi i/3)`, their slashes are `T:(x,y)->(x,omega y)` and `S:(x,y)->(-i/sqrt(3)) (x+2y,x-y)` by the same two- and three-term Poisson transform. The group has projective index 12 and cusp width 3. The valence formula bounds the order at infinity of a nonzero weight-27 holomorphic form by 27. The 28 monomials `x^(27-j) y^j` have distinct initial orders `j`, so triangular cancellation through `Q^27` proves that they span. Write `f=F(x,y)` for a homogeneous polynomial of degree 27, with coefficient `a_j` on the `j`th monomial.

Set `G_k=F(x+2y,omega^k(x-y))`, stripping off the common slash multiplier `(-i/sqrt(3))^27`. Its constant, from Poisson, is `3 sum_r n_r omega^(kr)`. Write `H_ij=[t^i](1+2t)^(27-j)(1-t)^j`. In the `x^(27-i)y^i` basis the `i=0..4` coefficients of `G_k` are its constant times `(1,0,0,-6,0)`: at `Q^3`, the contribution `27*6` from `x^27` is canceled against `3^3` from the third monomial. Inverting the three-character transform gives the **integer** equations

```
sum_(j=r mod 3) H_ij a_j = 3 n_r (1,0,0,-6,0)_i,  i=0,...,4; r=0,1,2.
```

The expansion `c + A8 Q^12 + O(Q^13)` first pins `a_0,...,a_12` by triangularity. For each `r` the five remaining `a_j` solve the five equations above (the exact replay checks the pivots). In particular, writing `G=G_0`, one gets

| Coefficient | Exact value |
|---|---|
| `[Q^13]f` | `354294 n_1` |
| `[Q^14]f` | `2499741 n_2` |
| `[Q^15]f` | `15707034 n_0 - 33 A8 - 398034 c` |
| `[Q^5]G` | `5235678 c + 27 A8 - 1194102 n_0 - 91854 n_1 - 367416 n_2` |
| `[Q^6]G` | `537345900 c - 567 A8 + 24596460 n_0 + 268272 n_1 + 4342653 n_2` |

One useful consequence for an individual nonzero coset of depth at least 8 (take `w_u=1`, all others zero) is that a class of type 1 has depth exactly `26/3` with 354,294 witnesses, and a class of type 2 has depth exactly `28/3` with 2,499,741 witnesses. In type 0, depth at least 10 forces depth exactly 10 with 15,707,034 witnesses. The type 2 count independently matches the known `Q+L` shell counted [from the elliptic curve](../uneven_layers55_20261122/LOCAL_OPTIMUM.md). These are not assertions that any additional such deep cosets exist.

## Applying positivity to all intermediate lattices

For `U=M/L ≤ A` take its indicator, so `c=1`, `n_0+n_1+n_2=3^s`. Orthogonality of characters makes `G=3*3^s sum_(v in U^perp) theta_v`; the sum is the theta series of `sqrt(3)M*` with variable `Q`. Both `Q^5` and `Q^6` coefficients are nonnegative. Their combination `[Q^6]G+21[Q^5]G` cancels `A8` and gives

```
479682 n_0 + 1660662 n_1 + 3373083 n_2 <= 647295138.
```

* At `mu=8`, a gain requires `s>=6`. If the quadratic restriction on `U` is nonzero, `n_0 <= (5/9)*3^s`: diagonalize, and use one nonzero coordinate if rank 1, otherwise two nonzero coordinates, in which every specified value has at most five preimages among nine. The displayed left side is then at least `1004562*3^s`, already too large at `s=6`. If `U` is entirely isotropic, the display excludes `s>=7`. For `s=6`, the two individual coefficients give `32,046,840 <= A8 <= 32,571,720` and, after dividing by `3^7`, the scaled dual's `Q^6+21Q^5` is 136,080. Thus its `Q^5` coefficient is at most 6,480.
* At `mu>8` every isotropic subspace of `U` lifts to an even integral lattice of minimum at least 10. An isotropic plane is impossible: for nine entirely isotropic classes and `A8=0`, the `Q^5` coefficient is negative. Any five-dimensional quadratic space over `F_3` has a totally isotropic plane, counting the radical, so `s<=4`. This follows from the ordinary nondegenerate Witt index (at least `floor((rank-1)/2)` for positive rank) plus the radical dimension. There are only three available minima `26/3,28/3,10` since `M ⊂ L*` and `L ⊂ M`.
* At `mu=26/3`, a gain requires `s>=4`, thus equality; without an isotropic plane a four-dimensional restriction must be nondegenerate nonsplit, with `(n_0,n_1,n_2)=(21,30,30)`. The `Q^5` coefficient is negative.
* At `mu=28/3`, a gain requires `s>=2`, and `[Q^13]f=0` forces `n_1=0`. Any quadratic form of rank at least two represents 1, as can be seen on two diagonal coordinates; rank 0 would make the minimum even. Hence the restriction has rank 1, value 2 on its nonzero component. Its radical has dimension at most one because there is no isotropic plane, so `s=2` and `(n_0,n_1,n_2)=(3,0,6)`. Again `Q^5` is negative.
* At `mu=10`, both `n_1,n_2` vanish. The isotropic plane exclusion then leaves just `s=1`; its scaled-dual norm-10 count is forced to 183,708.

Unlike the integral cases, none of the nonintegral exclusions assume a line/plane of new classes whose induced distances were sampled. The quadratic restriction is that of the entire glue code `U`.

## Nonlinear motifs at separation 10

For a finite set `S0⊂A` of `m` distinct translates define `w_u=m^(-1) #{(a,b) in S0^2: a-b=u}`. Then `c=1`, the weights sum to `m`, and the Fourier coefficients `hat w_v=m^(-1)|sum_(a in S0) exp(2 pi i b_A(a,v))|^2` are nonnegative. Poisson therefore makes `G=3 sum_v hat w_v theta_v` coefficientwise nonnegative even with arbitrary **nonadditive** `S0`. At squared separation 10 the primal coefficients above force `A8=0` and `n_1=n_2=0`, hence `[Q^5]G=5235678-1194102m`, giving `m<=4`. If `m=4`, any sum of four third roots of unity has squared modulus at least one: writing multiplicities `(a,b,c)`, the modulus squared is `a^2+b^2+c^2-ab-bc-ca`, a positive integer since `a+b+c=4`. So `hat w_v>=1/4` for every `v`. The norm-10 shell of `J` would give `[Q^5]G/3>=1347192/4`, but the exact value for `m=4` is 153090. This contradiction gives `m<=3`.

If two shifts survive, their nonzero difference `u` yields a depth-at-least-10 isotropic class (the individual formula gives depth exactly 10). Negation gives the same depth for `-u`, and all differences in `{0,u,-u}` are among `0,±u`, because `u` has order three. Thus either no such enlargement exists, or the three-phase bound is attained by a lattice. This conclusion does not carry over to arbitrary real offsets outside `L*` or to an unrestricted periodic packing.

## Exact replay

Run `python verify.py` (standard library). It enumerates the relevant 2D theta coefficients exactly, checks the cusp jet interpolation on five input columns, all rank/type cases (even allowing restrictions not embeddable in this particular discriminant), the exact density gates and the nonlinear four-phase exclusion. The underlying theta transformation and valence argument are stated above; the inherited 54D and dual minimum references are separate. No shortest-vector search, numerical optimization, or new density lower bound is used.
