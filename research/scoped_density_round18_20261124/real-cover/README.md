# A second half-coset shell improves real covering in ranks 48, 72 and 96

Let \(L\) be an extremal even unimodular lattice in dimension \(n=24k\),
\(k=2,3,4\), with squared minimum \(2k+2\). Put \(D=2k\). Its **squared**
covering radius \(R\) satisfies
\[
 R\le B(D):=\frac{(\sqrt{1+36D}-1)^2}{27}.
\]

| Rank | Previous all-real bound | New all-real bound |
|---:|---:|---:|
| 48 | \(14/3\) | \(4.5154374386\ldots\) |
| 72 | \(22/3\) | \(6.9828948250\ldots\) |
| 96, conditional on existence | \(10\) | \(256/27=9.4814814815\ldots\) |

This concerns the covering radius of this particular class of lattices,
not a new upper bound on unrestricted packing density. No equality is claimed.

## Two-shell identities

The [previous half-torsion theorem](https://github.com/alexlegeartis/SpherePackings/blob/6696a479a8a55e440ff8e65f10702702fe7d599c/research/scoped_density_obstructions_round15_20261123/48-72-quarter-torus/HALF_BOUND.md)
gives half-coset depth at most \(D\), in half-integer increments,
and an antipodal spherical 2-design at depth \(D\).
For a half-coset of depth \(d_0=D-1/2\), put \(d_1=d_0+1\).
Write \(N_i\) for the populations of the two shells \(z^2=d_i\), and
\(A_{j,i}\) for their harmonic sums in an arbitrary homogeneous harmonic
polynomial of degree \(j\). Then
\[
 A_{2,1}=-20A_{2,0},\qquad A_{4,1}=4A_{4,0}.
\]

| Rank | \(N_0\) | \(N_1\) |
|---:|---:|---:|
| 48 | 108 | 34928 |
| 72 | 6272 | 2165248 |
| 96, conditional | 380632 | 136063200 |

At infinity take the signed lattice theta series for the relevant
half-character, of weight \(12k+j\). Its first \(k+1\) ordinary
\(q\)-coefficients are \((1,0,\ldots,0)\) for the scalar and zero for
harmonics. As in the half-torsion proof, divide by the half-integral
character factor \(\theta_3^2\theta_4^2\) and factor
\([Y(X+Y)]^{D-1}\) from the invariant polynomial, where
\(X=\theta_2^4\), \(Y=\theta_4^4\).
The remaining degree is \(d=2k+1+j/2\), with basis
\((X+2Y)^{d-2i}X^{2i}\) for \(0\le i\le\lfloor d/2\rfloor\).
The transformation \(S\) swaps \(X,Y\) up to a common minus, replaces
the theta factor by \(\theta_3^2\theta_2^2\) up to a common minus, and
gives the shifted theta by Poisson summation. The scalar sign is
\((-1)^{6k}=+1\); in a given harmonic degree the common sign is irrelevant
to the coefficient ratio. The three exact coefficient matrices have
respective nullities 0, 1, 1. [verify.py](verify.py) checks the shifted
coefficients and ratios.

## Geometry of the two shells

Take \(x\) attaining \(R\), a nearest lattice point \(w\) to \(2x\),
\(u=2x-w\), and \(H=u^2\le R\). The half-coset \(w/2+L\) is antipodal.
If its depth is at most \(D-1\), antipodality gives
\(R\le D-1+R/4\), already stronger than \(B(D)\).
If its depth is \(D\), the 2-design gives a shell vector with
\(|z\cdot u|\ge\sqrt{H/12}\), so
\[
 R\le D+H/4-\sqrt{H/12}
 \le\max(D,D+R/4-\sqrt{R/12}).
\]
The last inequality is convexity in \(H\). Solving proves \(R\le B(D)\)
in this case, because \(B(D)>D\) and \(3B(D)/4+\sqrt{B(D)/12}=D\).

It remains to rule out depth \(d_0=D-1/2\) for any maximizing \(x\) with
\(R\ge r\), where \(r=4.515,\ 6.982,\ 9.481\) in the three ranks.
These are all below \(B(D)\).
Let \(\delta=d_0+H/4-R\). Antipodality and the definition of \(R\)
give \(\delta\ge0\), \(|z\cdot u|\le\delta\) on the first shell,
and \(|z\cdot u|\le1+\delta\) on the second shell.
Let \(M_{j,i}=\sum_{z^2=d_i}(z\cdot u)^j\). In particular,
\[
 M_{4,1}\le (1+\delta)^2 M_{2,1}.
\]
Use the zonal harmonic polynomials in \(z\),
\[
 (z\cdot u)^2-\frac{Hz^2}{n},\qquad
 (z\cdot u)^4-\frac{6Hz^2}{n+4}(z\cdot u)^2
       +\frac{3H^2(z^2)^2}{(n+2)(n+4)}.
\]
The two identities give
\[
 M_{2,1}=c_2H-20M_{2,0},\qquad
 M_{4,1}=c_4H^2-C H M_{2,0}+4M_{4,0},
\]
where
\[
 c_2=(20N_0d_0+N_1d_1)/n,\quad
 C=(120d_1+24d_0)/(n+4),\quad
 c_4=\frac{6d_1c_2}{n+4}+
     \frac{3(4N_0d_0^2-N_1d_1^2)}{(n+2)(n+4)}.
\]
Thus the difference \(M_{4,1}-(1+\delta)^2M_{2,1}\),
divided by \(H>0\), is at least
\[
 c_4H-CN_0\delta^2-c_2(1+\delta)^2.
\]
Here we only required \(M_{2,0}\le N_0\delta^2\) and dropped
nonnegative terms.

Since \(R\ge\max(r,H)\),
\(0\le\delta\le\min(d_0+H/4-r,d_0-3H/4)\),
so \(H\) is in \([4(r-d_0),4d_0/3]\). Substitute that upper bound on
\(\delta\) into the last displayed lower bound. On
\([4(r-d_0),r]\) the result is concave in \(H\), and on
\([r,4d_0/3]\) it is increasing because \(c_4>0\) and the substituted
nonnegative \(\delta\) decreases. Hence only the values at
\(H=4(r-d_0)\) and \(H=r\) matter. They are both positive, as verified
over the rationals in [verify.py](verify.py), contradicting the required
upper bound on \(M_{4,1}\). This proves the claimed full real-torus bound.

## Consequence for arbitrary real shifts in full-layer stacks

Consider any periodic packing with complete horizontal translates of
\(L\), in \(n+s\) dimensions, at squared separation \(\mu>B(D)\).
There is at most one phase at each vertical site and distinct vertical
sites are at squared distance at least \(\mu-B(D)\). If \(c_s\) bounds
the vertical center density (general periodic or lattice, as appropriate),
the stack's center density is at most
\[
 \frac{c_s}{2^n}\mu^{(n+s)/2}(\mu-B(D))^{-s/2}.
\]
This follows by bounding the vertical point density by
\(c_s(2/\sqrt{\mu-B(D)})^s\); each complete horizontal layer has
covolume one. It applies to all real phases, not only specified torsion
classes. At the parent minimum \(\mu=D+2\), it improves the older
real-phase upper bound for the 5-extra-dimensional 48D stack
by about 24% and the 4-extra-dimensional 72D stack by about 57%.
For the hypothetical 96D parent, the older all-real estimate gave
no finite ceiling at \(\mu=10\); the new one does.
This does not close the gap to the known half-torsion constructions.

Replay: python verify.py (Python 3 and SymPy; exact arithmetic only).

## Attainment restrictions

In an integral lattice there is a rational Voronoi vertex realizing \(R\):
the active independent equations are \(x\cdot v=v^2/2\) in a lattice basis.
In particular \(R\) is rational, so the inequalities in ranks 48 and 72
are strict for each individual parent.

If the conditional rank-96 endpoint \(256/27\) is attained, each doubled
iterate is also maximizing, and each chosen half shell is at depth 8.
Indeed the shallower branches have already been excluded, and convexity
is strict for \(0<H<R\) in the deepest branch. The sharpness of the
second-moment step then forces every contact to have
\(|z\cdot u|=8/9\). Orient the \(8,222,400/2\) antipodal pairs by
positive projection; the sum \(W\) of the oriented representatives lies
in \(L\), as there are an even number from one half-coset.
The frame identity gives
\[
 W=\frac{8,222,400}{24(8/9)}u=385425u.
\]
Thus each doubled maximizing phase has exact order dividing
\(385425=3^3\cdot 5^2\cdot571\). It must be divisible by 9, since
\(256m^2/27\) must be an even integer for any such order \(m\).
This is a necessary condition, not a rank-96 attainment or existence
claim.
