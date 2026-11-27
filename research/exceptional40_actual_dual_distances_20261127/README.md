# Actual-distance optimum of all dual-coset motifs on the exceptional 40D carrier

**No density record.** Start with any extremal even unimodular 48-dimensional lattice containing a primitive normal `K=sqrt(3)E8`; let `L` be the rank-40 orthogonal section. Among **all unions of full translates of L using phases in L*/L**, with arbitrary nonlinear phase choices and any packing radius, the center density is at most `3^16/2^20`. Equality comes only from one phase. This is a statement about the **actual horizontal coset minima**, unlike the [standard normal-diameter bound](../exceptional_torus_40_42_20261126/). It does not cover arbitrary real phases, except for the local deep-phase statement below, or other carriers.

The **unrestricted** named comparison is Dorofeev–Sun–Wang, [*Optimal Extensions of Cross-Sections*, Table 1, §3 and Appendix A.1](https://arxiv.org/html/2607.20359v4); their 40D packing is a lattice and is the singleton case here. The [live project table](../../LOWER_BOUNDS_1_96.md) gives about `1.47983399313426e-7` as an occupied fraction. The conversion is multiplication by `pi^20/20!`. We use the exact center density, not that rounded number, in all comparisons. The argument requires no assumptions on which of the archived embeddings is used and no ambient shortest-shell enumeration.

## A low-degree model for the actual coset theta functions

The unimodular primitive gluing makes `A=L*/L` anti-isometric to `K*/K` for the discriminant quadratic form `q_L(u)=||u||²/2 mod Z`. In particular `det L=3^8`, `A≅F3^8`, and `J=sqrt(3)L*` is even. The class minima of `K*/K` are at most `3`: `sqrt(3)E8` has squared covering radius 3. For any `x in L*` outside `L`, choose a matching normal representative `y` of squared norm ≤3. The vector `(x,y)` is a nonzero vector of the ambient lattice, so `||x||²≥6−3=3`. As `3||x||²` is even, `||x||²≥10/3`. Nonzero `x in L` already have squared norm ≥6. Thus **J has squared minimum at least 10**. This is the source of the cusp vanishing in the calculation below.

Use `q=exp(2 pi i z)`, `Q=q^(1/3)`, and `theta_u(z)=Σ_(x in L+u) q^(||x||²/2)`. The coset theta-vector is a weight-20 holomorphic form for the level-three Weil representation, trivial on `Gamma(3)` ([Scheithauer, Theorem 4.7 and the principal-congruence triviality used there](https://doi.org/10.1093/imrn/rnn166)). Its `T`-action multiplies `theta_u` by `exp(2 pi i q_L(u))`. Under weight-20 slash by `S`, Poisson gives `(-i)^20/81` times the sum of `theta_v` weighted by `exp(-2 pi i b_L(u,v))`. Reversing the Fourier sign would not change any constant or coefficient used here. In particular the coefficients of `Q^1,...,Q^4` vanish at all four cusps represented by the identity and `T^k S` for `k=0,1,2`. This holds for every linear combination of the `theta_u`.

Let `x=theta_(A2)` and `y=theta_(A2+(1/3,1/3))` in a basis of Gram `[[2,1],[1,2]]`. They are weight-one forms on `Gamma(3)`, with `x=1+6Q^3+O(Q^9)`, `y=3Q+3Q^4+O(Q^7)`. For `omega=exp(2 pi i/3)`, their slashes are `T:(x,y)→(x,omega y)` and `S:(x,y)→(-i/sqrt(3))(x+2y,x−y)`. These follow from the three-class Poisson transform, since the two nonzero cosets are negatives. In any positive integral weight `w` used here, the monomials `x^(w−j)y^j`, `0≤j≤w`, have distinct initial orders, and span the holomorphic space: a form with vanishing order ≥`w+1` at infinity would contradict the valence bound `w`, since the projective index and the cusp width are 12 and 3, respectively.

For any real weights `w_u` put `f=Σ w_u theta_u=F(x,y)`, `c=w_0`, `A5=[Q^5]f`, `n_r=Σ_(q_L(u)=r/3) w_u`. Write `a_j` for the coefficient of `x^(20−j)y^j` in `F`. The coefficients through `Q^5` at infinity pin `a_0,...,a_5`. At the other cusps strip off the common factor `(-i/sqrt(3))^20`, so the expansion there is `G_k=F(x+2y,omega^k(x−y))`. Its constant term is `729 Σ_r n_r omega^(kr)`, the ratio between the two Poisson constants. Its first five coefficients in the `x^(20−i)y^i` basis are that constant times `(1,0,0,−40/9,0)`: the `Q^3` term is canceled since `20·6−(40/9)·27=0`. The finite character transform thus gives the equations

```
sum_(j=r mod 3) H_ij a_j = 729 n_r (1,0,0,−40/9,0)_i,
H_ij=[t^i](1+2t)^(20−j)(1−t)^j,      i=0,...,4; r=0,1,2.
```

The five unknown coefficients in each residue class solve five exact linear equations. The accompanying replay uses integer `A2` theta enumeration and rational elimination. In the resulting identity `G=G_0`, the relevant coefficient values are:

| Coefficient | Exact expression |
|---|---|
| `[Q^6]f` | `495(n0−c)` |
| `[Q^7]f` | `9180 n1` |
| `[Q^8]f` | `116235 n2−40 A5` |
| `[Q^9]f` | `1084590 n0+7410 c` |
| `[Q^5]G` | `22320522 c +59049 A5 −275562 n0 +551124 n1 −1102248 n2` |

For `f=theta_0`, `c=n0=1` and `A5=n1=n2=0`. It follows that `L` has exactly **1,092,000** vectors of squared norm 6, verifying its minimum as well as the prior contact count. Also `G=729 Σ_u theta_u`, so `[Q^5]G/729=30,240` is the number of norm-10 vectors of `J`. Thus `J` really attains 10 and at least one nonzero horizontal coset has squared depth `10/3`.

For a single nonzero phase, take `c=0` and the appropriate `n_r=1` with the other two zero. The allowed `Q` exponents in that coset are `r mod 3`. The actual squared minima are therefore:

| `q_L(u)` | Squared depth | Number of nearest vectors when specified |
|---|---:|---:|
| `0`, `u≠0` | `4` | `495` |
| `1/3` | `14/3` | `9180` |
| `2/3`, `A5>0` | `10/3` | `A5` |
| `2/3`, `A5=0` | `16/3` | `116235` |

This does not assert equality with the conservative normal formula separately for every long class. In particular the distribution of the 30,240 points across individual type-2 classes is not needed. There are deep phases: a matching normal class of squared depth `2/3` forces the horizontal depth at least `16/3`, since the ambient norm is at least 6.

## All nonlinear motifs and all radii

Let `m>1` distinct shifts form the period cell. Its volume per point is `81/m` and the full squared minimum `mu` is one of the four depths above. Compared to the singleton, the center-density ratio is `m (mu/6)^20`.

* **`mu=16/3`.** Every pair difference has `q_L=2/3`, equivalently the matching normal finite quadratic form `3q_K` takes the value 1 in `F3`. Translate a site to zero. On any nine other site vectors, the Gram of the polar form over `F3` would be `I_9+J_9`: diagonal value 2 and off-diagonal value 1. This matrix is invertible since `1+9≠0 mod 3`, impossible in dimension eight. Thus `m≤9`.
* **`mu=14/3`.** Nonzero isotropic differences cannot occur. The normal discriminant space `K*/K≅E8/3E8` has an isotropic subgroup of size 81. For explicit coordinates, use doubled E8 coordinates reduced modulo 3: they give all `F3^8` classes since the coordinate denominator is a power of 2. The polar form of `3q_K` is ordinary dot product mod 3, since `4=1 mod 3`. In each disjoint four-coordinate block choose the pair `(1,1,1,0)`, `(1,2,0,1)`, which are self-orthogonal and mutually orthogonal over `F3`; the four rows span the required subgroup. Each of its 81 cosets can contain at most one site. Hence `m≤81`.
* **`mu=4`.** Some nonzero difference `v` has actual depth `10/3`, by the total scaled-dual count above. All pairs in a coset of `<v>` differ by `±v` and are incompatible. There are `3^8/3=2187` such cosets, so `m≤2187`.
* **`mu=10/3`.** The full group has size 6561, so `m≤6561`.

Each of these four caps paired with its radius gives `m (mu/6)^20<1` by exact rational powers. Any smaller chosen radius only lowers density. This handles arbitrary finite motifs rather than only subgroups, and in particular excludes every improving integral or nonintegral intermediate lattice `L ⊂ M ⊂ L*`.

## No small off-grid slip of the deepest phases

For a homogeneous harmonic polynomial `P` of degree `d=1` or 2, the coset theta with inserted factor `P(x)` has weight `20+d` and the same finite `T` eigenvalue. Harmonic Poisson expresses its other cusps as linear combinations of the weighted coset series, with a nonzero common scale. At **all four cusps** its expansion has order at least 5: the zero vector contributes 0 because `d>0`, and all other norms in `L*` have `Q` exponent at least 5. In the degree-`20+d` homogeneous `x,y` polynomial this forces divisibility by `y^5 (x^3−y^3)^5`. Indeed at infinity the factor is `y^5`; at the three others apply the explicit slashes above, whose vanishing at the transformed cusp points gives the three factors `x−omega^k y`, each to the fifth power. The product's `T` eigenvalue is `omega^2`.

For a phase with `q_L=2/3` and depth `16/3`, the residual polynomial has degree one or two and `T` eigenvalue 1. It must be a multiple of `x^d`. If the multiple were nonzero the weighted theta would start at `Q^5`, but this coset has no such vectors. Thus the entire harmonic theta is zero for degrees 1 and 2. Its set of 116,235 nearest displacements has mean zero and second moment `(16/3) I_40/40`, i.e. it is a spherical **2-design** at that radius.

This makes every such phase a strict local maximum of squared distance to `L`. More explicitly, for any displacement `s` the nearest displacements `v` have `mean(v·s)=0`, `mean((v·s)^2)=(16/3)||s||²/40` and `v·s≤sqrt(16/3)||s||`. For values in `[-a,b]` with zero mean, their mean square is at most `ab`; thus some `v·s≤−sqrt(16/3)||s||/40`. Therefore for `0<||s||<sqrt(16/3)/20`, the depth at the shifted phase is strictly smaller than `16/3`. This is local only: in particular the existence of a real off-palette phase of depth `>6·2^(−1/20)` (which would allow a two-phase improvement) is not addressed.

Run `python verify.py`. The standard-library replay verifies all coefficient identities by checking the five linear input columns, the finite isotropic subgroup, cap arithmetic and the harmonic residue check. There is no high-dimensional search or unverified minimum claim in the proof.
