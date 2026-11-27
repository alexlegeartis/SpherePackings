# The same 9D result for squared minimum between 29/25 and 2

Keep the [same horizontal E7 and arbitrary colored planar Bravais
normal](README.md). Fix the horizontal minimum at 2, but let the full
packing's squared minimum be \(\mu\in[29/25,2]\). Its center density is
still at most \(1/(16\sqrt2)\); it is strictly smaller when \(\mu<2\).
Lower \(\mu\), partial fibers and nonlattice normal projections
are not covered.

The main proof's root descent works on **every** dilate \(\alpha R\):
if \(x\cdot r>1\), compare with \(z-(r\cdot z)r\) as before, but now
the polar bound is \(x\cdot z\le\alpha\). After the descent,
\(x\in\alpha R\cap Q\) gives \(\|x\|^2\le\alpha\). In particular, a
horizontal phase represented in \(\alpha R\) has torus square at
most \(\alpha\).

Let \(a\) be the squared length of a shortest nonzero normal vector
\(A\). If \(a<\mu-1\), look at the consecutive fibers at \(0,A,2A\).
Both horizontal differences have torus square at least \(\mu-a>1\).
By the root-polar bound (1), both lie modulo \(E_7\) in
\(h+(a+3/2-\mu)R\). As \(2h\in E_7\), the phase difference from
the first to the third fiber is represented in
\(2(a+3/2-\mu)R\). Its torus square is thus at most
\(2(a+3/2-\mu)\), and
\[
4a+2(a+3/2-\mu)\ge\mu,\qquad a\ge(\mu-1)/2.
\]
If \(a\ge\mu-1\) initially, the latter bound already holds.

Write \(d=2-\mu\), and multiply the squared metric of the normal plane by
\(\lambda=1+d/a\), holding the horizontal geometry and every phase
fixed. Every inter-fiber squared distance increases by at least \(d\);
intrafiber distances already have minimum 2. The planar covolume is
multiplied by \(\lambda\), so the main theorem says the original
density divided by the incumbent is at most
\[
\lambda(1-d/2)^{9/2}\le
f(d)=\frac{1+d}{1-d}(1-d/2)^{9/2}.
\]
This is at most one for \(0\le d\le21/25\): the logarithmic derivative
has the sign of \(9d^2-4d-1\), hence there is at most one turning point,
a minimum, on the interval. At \(d=0\), \(f(d)=1\); at \(d=21/25\),
\(f(d)^2=(23/2)^2(29/50)^9<1\) by direct integer comparison.
