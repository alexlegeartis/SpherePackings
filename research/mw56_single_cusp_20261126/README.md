# MW56: full one-cusp fibers cannot improve 57D or 58D

**Exact scope.** Let $L=MW_{56}$ be the even Mordell–Weil lattice of
$y^2=x^3+t^{30}+1$ over $\mathbb F_{29^2}(t)$. Fix any one finite cusp.
In each affine horizontal plane, use an entire fiber of the additive reduction
$L\to\mathbb F_{29^2}$ at that cusp. Any collection of such fibers stacked in
one or two perpendicular dimensions has at most **54.8%** or **35.3%** of the
respective existing 57D or 58D density, at any packing radius. Normal positions
and color assignments need not form a lattice or be periodic.
This does not concern arbitrary translations of $L$, arbitrary proper
subsets of reduction fibers, or using several cusps at once. No record changed.

## 1. The exact distance table of cusp fibers

Use field elements $a+bs$, $s^2=2$, encoded as $a+29b$ with
$0\leq a,b<29$. The curve has smooth infinity, thirty irreducible type-II
fibers and $\chi=5$. By the usual height formula the lattice norm is
$10+2(P\cdot O)$ for nonzero sections. A point represented by polynomials of
degrees at most 10 and 15 is disjoint from $O$ and has norm 10.
The parent rank 56 and determinant $29^8$ are the $p=29$ instance of
[Schütt–Shioda, *Mordell–Weil Lattices*, Theorem 13.35](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf);
the [previous exact construction](../packing_binary57_58_20261117/) uses
this identical parent. Conversely, norm 10 means $P\cdot O=0$, so its
coordinates are polynomials of those bounded degrees.

At a cusp $t_0^{30}=-1$, reduction to the smooth locus of
$y^2=x^3$ identifies that group with the additive field, via
$\rho(P)=x(t_0)/y(t_0)$ on nonidentity reduction and $\rho(P)=0$ on identity
reduction.
(In the smooth affine locus $x=1/\rho^2,\ y=1/\rho^3$; the chord-and-tangent
sum is additive in $\rho$.) No norm-10 point reduces to zero: its
coordinates are finite, and $x,y$ cannot simultaneously vanish where
$t^{30}+1$ has a simple zero. Thus $\ker\rho$ has minimum at least 12.

The finite-field unitary group acts by linear substitution in the homogeneous
equation $Y^2-X^3=T^{30}+S^{30}$ (degrees 15, 10, 1, 1).
It is transitive on the isotropic lines, and hence on the cusps.
For an isotropic line spanned by $v=(t_0,1)$, choose a second isotropic
vector $w$ with Hermitian pairing 1: start with any vector having pairing
1 and subtract half its self-pairing times $v$. In this basis
$\operatorname{diag}(\lambda,\bar\lambda^{-1})$ is a unitary map for
any $\lambda\in\mathbb F_{29^2}^{\times}$. Thus its stabilizer realizes
every eigenvalue $\lambda$ on $v$. In affine coordinates a unitary
matrix acts by $t\mapsto(at+b)/(ct+d)$; section coordinates get multiplied
by $(ct+d)^{10}$ and $(ct+d)^{15}$. At the fixed cusp, $ct_0+d=\lambda$,
so the cusp color is multiplied by $\lambda^{-5}$. It follows that
five point witnesses with different values of $\rho^{168}$ suffice
to cover all $29^2-1=840$ nonzero elements by this action.

Here are such witnesses at $t_0=12$, written as
$x=at^{10}+bt^5+c,\ y=ut^{15}+dt^{10}+et^5+f$.

| $a,b,c;u,d,e,f$ | $\rho$ | $\rho^{168}$ |
|---|---:|---:|
| 0, 3, 0; 1, 0, 0, 28 | 24 | 1 |
| 2, 416, 483; 3, 832, 579, 290 | 231 | 84 |
| 2, 425, 396; 3, 821, 318, 290 | 41 | 365 |
| 2, 454, 483; 3, 38, 579, 551 | 824 | 510 |
| 2, 29, 20; 3, 58, 2, 290 | 337 | 809 |

Therefore **every nonzero color coset has minimum exactly 10**, and $\rho$
is onto, giving covolume $29^4\cdot841=29^6$ for its kernel. The ancestry
of the witnesses is the rational elliptic surface $Y^2=X^3+u^6+1$ under
$u=t^5$; in contrast, the simple sections $x=-(at+b)^{10},\ y=(ct+d)^{15}$
alone give just one of the five multiplicative classes at a cusp.

For the other end of the table, the certificate contains two
polynomial points of the same color $23$ at $t_0=12$. Their
coordinate-difference gcd is $t+17$ and they do not meet at infinity.
Their affine meeting is in the smooth locus of the cusp fiber, so
$P\cdot Q=1$, including multiplicities. The pairing is $5-1=4$, hence
$P-Q\in\ker\rho$ has squared norm $10+10-2\cdot4=12$.
By transitivity this shows **the zero-color lattice has exact minimum 12**
at each cusp as well.

## 2. No-record theorem, independent of the normal pattern

Write $Q=841$, $k=1,2$, $d=56+k$, and let $R$ be the desired
squared center separation; necessarily $R\leq12$ because each full fiber
contains a pair at separation $\sqrt {12}$.
The set of positions of each given color has normal minimum separation
$\sqrt R$, since two identical fibers have a horizontal pair with zero
difference. For $R>10$, distinct colors cannot share a normal position
and the union of all normal positions has separation at least
$\sqrt{R-10}$, since a difference of distinct fibers has a horizontal
vector of squared norm 10. These are necessary constraints, not assumptions
about additivity of the normal placement.

Thus the upper number density per normal volume is bounded by
$$
 A_k(R)=c_k\min\{Q R^{-k/2},(R-10)^{-k/2}\},\qquad
 c_1=1,\quad c_2=2/\sqrt3,
$$
with only the first term used for $R\leq10$. These constants are the
1D interval and 2D triangular packing bounds. Horizontal number density
is $29^{-6}$ per fiber, uniformly in its coset; taking larger and larger
product boxes gives the same bound for upper (not necessarily periodic)
packing density. The center density is at most $(R/4)^{d/2} A_k(R)/29^6$.

Put $R_*=10/(1-Q^{-2/k})$. Below $R_*$ the bound is increasing.
Above $R_*$ it is proportional to $R^{d/2}(R-10)^{-k/2}$, whose
logarithmic derivative has the sign of $56R-10d$. Hence it has no
interior maximum and only $R_*$ and $R=12$ need comparison.
The existing exact 57D/58D comparators have center densities
$(5/2)^{57/2}/(\sqrt3\,29^4)$ and
$(5/2)^{29}/(\sqrt{27/4}\,29^4)$, from the
[independently audited half-coset extensions](../packing_binary57_58_20261117/).
The corresponding ratios, unchanged for ordinary density, are:

| Dimension | At $R_*$ | At $R=12$ | Uniform upper ratio |
|---|---:|---:|---:|
| 57 | 0.547744241314 | 0.262975413151 | $<0.548$ |
| 58 | 0.310162384857 | 0.352818540130 | $<0.353$ |

Replay the literal polynomial identities, kernel witness, and all exact
ratio comparisons using python3 cusp_color_certificate.py (stdlib only;
do not disable assertions).
