# A sharp 23D layer ceiling for the \(\sqrt2D_2\) Leech section

**Result.** Let \(K\) be the rank-22 section of the Leech lattice perpendicular to the primitive \(\sqrt2D_2\) supported in its first two standard Golay coordinates. Stack whole translates of \(K\) at arbitrary real heights in one new direction, restricting the horizontal phases to the 16 classes of \(K^*/K\). The sharp center-density ceiling is **\(1/2\)**, even with nonuniform heights, nonlinear ordering, and multiple phases at one height. A period-four lattice stack attains it. This is below the **external** 23D benchmark \(23^{23/2}/(2^{34}3^{12})\) from [Chen, Hu, Li, Wang and Wu (2025), §3](https://arxiv.org/abs/2505.02394). It is not a global 23D bound.

The structural step is that a dense run of layers cannot keep using locally cheapest moves: the color choices would close a four-cycle and place two copies of the same phase too close.

## Exact phase geometry

Use raw integer coordinates, divided by \(\sqrt8\), for the Leech lattice \(\Lambda\). The first two coordinates of its orthogonal projection have the same parity, and the supported section in their plane \(U\) is

\[
N=\tfrac1{\sqrt8}\langle(4,4),(4,-4)\rangle \cong\sqrt2D_2.
\]

Indeed, a vector supported there must be even (the other 22 coordinates vanish); its residue mod 4 cannot use a nonzero Golay support of size at most two; and its sum is divisible by eight. This gives exactly the displayed section. By primitivity and unimodularity, \(K=\Lambda\cap U^\perp\) has covolume \(4\), and normal classes \(N^*/N\) correspond bijectively to horizontal classes \(K^*/K\) by Leech gluing. In raw coordinates \(N^*=\{(p,q):p\equiv q\bmod2\}/\sqrt8\).

For each nonzero class \(s\), let \(\gamma(s)\) be its shortest normal squared norm. Each representative needed below occurs as the first two raw coordinates of a norm-4 Leech vector. This is checked from the three familiar Golay shapes: two entries \(\pm4\); eight entries \(\pm2\) on an octad with even number of minus signs; or the Golay sign changes of a single \(-3\) among 23 entries \(1\). No global closest-vector search is needed. For any lift of a horizontal vector in this class we can subtract a vector of \(N\) to obtain the shortest normal component; the ambient norm is at least 4 and the shell witnesses attain equality. Thus the horizontal squared depth is exactly \(d(s)=4-\gamma(s)\).

| Shortest raw normal shapes, up to signs/order | Number of classes | \(\gamma\) | \(d\) |
| --- | ---: | ---: | ---: |
| zero | 1 | 0 | 0 |
| \((1,1)\) | 4 | \(1/4\) | \(15/4\) |
| \((2,0)\) | 4 | \(1/2\) | \(7/2\) |
| \((2,2)\) | 2 | \(1\) | 3 |
| \((3,1)\) | 4 | \(5/4\) | \(11/4\) |
| \((4,0)\) | 1 | \(2\) | 2 |

The zero depth is the minimum horizontal separation between two copies of the same phase at distinct heights. Within each individual sheet the lattice minimum is 4. We use the actual minimum \(\mu\) of the stack, so \(0<\mu\le4\). If its mean number of sheets per unit vertical length is \(R\) (counting different phases at one height), then its center density is \((\mu/4)^{23/2}R/4\).

## Phase anticode bounds

Let \(T\) be a subset of the 16 classes whose pairwise horizontal depths are bounded by \(d_0\). All heights carrying classes in any translate of \(T\) are separated by at least \(\sqrt{\mu-d_0}\) if \(\mu>d_0\). Average over the 16 translates of \(T\): each phase is counted \(|T|\) times. Consequently \(R\le 16/(|T|\sqrt{\mu-d_0})\).

Use the following explicit anticodes:

1. A singleton has \(d_0=0\).
2. The four classes \((\pm2,0),(0,\pm2)\) have \(d_0=3\).
3. The eight classes with even raw coordinates have \(d_0=7/2\).
4. The entire group has \(d_0=15/4\).

These bound the center density by \(1/2\) over \((0,7/2]\), \([7/2,37/10]\), \([37/10,379/100]\), and \([31/8,4]\), respectively. Each bound is a constant times \((\mu/4)^{23/2}/\sqrt{\mu-d_0}\). Its logarithmic derivative has the sign of \(22\mu-23d_0\), so the endpoint comparisons are sufficient and are rational after squaring.

## Four gaps cannot all use the cheapest route

In the remaining interval \([379/100,31/8]\), all heights are distinct. Set

\[
a=\sqrt{\mu-15/4},\quad b=\sqrt{\mu-7/2},\quad
c=\sqrt{\mu-3},\quad e=\sqrt{\mu-11/4}.
\]

Each of four successive gaps is at least \(a\). If some phase increment is not one of the four signs \((\pm1,\pm1)\), its gap is at least \(b\). Otherwise, if two adjacent increments do not toggle exactly one sign, they are equal or opposite, so their sum has depth at most 3 and the two intervening gaps sum to at least \(c\). If each neighboring pair does toggle exactly one sign but a triple toggles the same sign twice, its sum has normal class of type \((3,1)\) and the corresponding gaps sum to at least \(e\). If none of these happens, the four sign pairs run once around their square and sum to zero, so their gaps sum to at least \(\sqrt\mu\). In all cases the four-gap total is bounded below by

\[
M(\mu)=\min\{b+3a,\ c+2a,\ e+a,\ \sqrt\mu\}.
\]

Averaging overlapping windows gives \(R\le4/M(\mu)\), so the density is at most \((\mu/4)^{23/2}/M(\mu)\). To certify this at most \(1/2\) everywhere in the interval, the verifier partitions it into intervals of width at most \(1/200\), uses the left endpoints in each increasing denominator and the right endpoints for the numerator, and compares rational outward bounds after squaring. This arithmetic is independent of any packing search.

## Attainment and remaining route

At heights \(j/2\), take the phase \(j(1,1)\). Gaps of one and two indexed layers have total squared separation \(15/4+1/4=4\) and \(3+1=4\); three steps are farther, and four or more have vertical squared length at least 4. The density is \(1/2\) and additive gluing makes it a lattice. This does not claim new priority.

The exact restriction is to full \(K\) sheets and its 16 discriminant phases. Arbitrary continuous horizontal phases, modified horizontal stock, or a different primitive cut are not excluded. We also tested a distinct, explicitly structured 3-translate geometry on the norm-6 Leech section: a 96-by-96 rational grid of the two-dimensional torus in which coordinates are constant on a Golay octad and its complement, apart from the distinguished coordinate. No grid triple approached the \( (23/6)(2/3)^{2/23}\approx3.7005\) squared-depth threshold; the best triple among phases deep enough from the origin had third depth \(10/3\), while the individual maximum observed was \(23/6\). That finite probe is not a continuous upper bound.

Run the companion Python verifier for the exact class census, shell lifts, sample attainment, and the rational density comparisons. See also the independent [shortened-Golay proof](GOLAY.md) and its exact finite ingredients, which close another two-translate carrier without assuming the shift is rational.
