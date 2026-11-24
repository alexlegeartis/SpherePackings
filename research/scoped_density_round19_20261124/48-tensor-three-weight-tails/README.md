# Three exceptional E₈ slots: the large-weight negative tail

In the setup of the three-exception theorem, if \(a<0\) and \(c\ge129\),
then the same universal scoped center-density cap 16,800 holds.
This argument is independent of the motif and does not use any Fourier
computation. It does not address \(a<0,1<c<129\).

Let \(D=ab+ac+bc+4abc\), \(H=11/4000\), and \(F=b+a(1+4b)\).
Positive definiteness of a fixed-\(c\) slice gives \(F>0\);
also \(z=1+4a>0\) by using a period vector supported on the
negative slot and the four ordinary slots. In particular \(|a|<1/4\).
As before, the rank-three Hermite bound closes \(D\le H\).
Assume \(D>H\) below.

Use the period sublattice whose integer coefficients along a single E₈
root are \((-(u+4r),u,0,r,r,r,r)\). The two-dimensional scalar form
before multiplication by 4 has Gram determinant \(4F\). The exact
planar Hermite constant gives \(g\le4\sqrt{16F/3}\). The period vector
with \(u=0,r=1\) gives \(g\le16z\). We also have \(g\le4(1+a)<4\).

If \(D\ge cF/2\), the trivial phase bound gives
\[
 \rho\le 2^{48}\frac{\min(2^{-24},(4F/3)^{12})}{(cF/2)^4}.
\]
The bound increases up to \(F=3/16\) and decreases thereafter. Thus
\(\rho\le2^{24}(32/(3c))^4<16800\).

If \(D<cF/2\), from \(D=cF+ab\) obtain \(F<2|ab|/c<b/(2c)\) and
\(z=(F-a)/b<1/(2c)+1/(4b)\). Hence
\(g\le4\sqrt{8b/(3c)}\) and \(g\le8/c+4/b\).
If \(b\le19/5\) take the first estimate; otherwise take the second.
For \(c\ge129\), the square of this bound is at most
\[
 Q=\max\!\left(\frac{128(19/5)}{3\cdot129},
            \left(\frac8{129}+\frac4{19/5}\right)^2\right).
\]
Using \(D>H\) and at most \(2^{48}\) phases gives
\(2^{48}(Q/64)^{12}/H^4<16800\), as verified by exact rational replay.

## A second tail: both heavy weights at least 7

Still with \(a<0\), if \(b\ge7\) (and therefore \(c\ge7\)), then
the same cap holds even for \(c<129\). In fact
\[
 D/(bc)=z+a(1/b+1/c)\ge (15z-1)/14,
\]
so \(D\ge 7(15z-1)/2\) for \(z\ge1/15\).
Let \(z_0=(1+2H/7)/15\). For \(z\le z_0\), use \(g\le16z_0\) and \(D>H\).
For \(z_0\le z\le1/4\) use \(g\le16z\) and \(D\ge 7(15z-1)/2\).
The upper bound in this interval is proportional to
\(z^{24}/(15z-1)^4\). Its only possible interior turning point is
a minimum (the derivative has the sign of \(300z-24\)).
For \(z\ge1/4\) use \(g\le4\) and the increasing determinant lower
bound. Thus only \(z=z_0\) and \(z=1/4\) need to be tested.
Exact comparisons are in the same verifier.

Together these two tails confine any outstanding case to
\(a<0,\ 1<c<129,\ b<7\).

## Moderate heavy weight and no near-cancellation

There is a stronger cutoff in the first determinant branch: if
\(c\ge4\) and \(D\ge cF/2\), the cap follows independently of
the size of \(b>0\). Fix the \(c\) slot in each phase fibre.
The remaining six zero-sum \(E_8\) slots are a periodic
40-dimensional packing with period determinant \(F\) on the
five independent scalar slots. Its center density is
\(N_{\rm fibre}(g/8)^{20}/F^4\), hence bounded by
\[
 B=\frac{42\cdot 20!}{2^{21}\,3^{20}}<16800.
\]
For completeness, for any packing of minimum distance 2 the
gauge \(f(x)=(1-\|x\|^2/2)_+\) has total contribution at
each location at most 1. For \(m\) centers where the
gauge is positive, summing their pairwise squared distances
gives \(\sum\|x_i\|^2\ge 2(m-1)\), proving the assertion.
Integration gives the center-density cap
\(42/(2^{21}V_{40})\le B\) by \(\pi>3\).
This argument is invariant under scaling.

With \(256\) possible fixed labels the 48D center density is bounded by
\(256B(gF/(8D))^4\le 256B/c^4\le B\), using \(g<4\)
and \(D\ge cF/2\).
