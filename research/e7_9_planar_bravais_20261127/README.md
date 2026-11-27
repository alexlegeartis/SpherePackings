# A sharp 9D bound for complete E7 fibers on a planar lattice

Let horizontal \(E_7\) have minimum squared length 2. Suppose the projection
of the fiber positions onto the orthogonal two-plane is a full Bravais
lattice, with **one complete translate of \(E_7\) at each site**. Its
translation may be arbitrary, continuous and nonperiodic as a function of
the site. If the full packing has minimum squared distance 2, its center
density is at most the classical \(\Lambda_9\) value \(1/(16\sqrt2)\).
This is a restricted result, not a proof of 9D optimality. It does not
cover nonlattice normal site sets, partial \(E_7\) fibers, or smaller
squared minimum for the same horizontal \(E_7\).

## 1. Root-polar geometry

Use the lattice \(L=E_7\) in the zero-sum subspace of \(\mathbb R^8\):
the coordinates are either all integral or all half-integral. The two
classes differ by the half-vector with four plus and four minus signs,
so \(L\) has covolume \(\sqrt8/2=\sqrt2\), and its minimum is 2.
Its roots (squared norm 2) are \(e_i-e_j\) and the half-vectors with
four \(+1/2\) and four \(-1/2\). Put \(Q=\{x:x\cdot r\le1\ \forall r
\text{ a root}\}\) and \(R=\operatorname{conv}\{r:r\text{ a root}\}\).
The centrally symmetric bounded polytopes \(Q,R\) are polar. Every
nearest representative of a torus point modulo \(L\) belongs to \(Q\),
by comparison with the roots. No converse or Voronoi-relevance
theorem is needed.

For sorted coordinates, the root inequalities become
\(x_1-x_8\le1\) and \(x_1+\cdots+x_4\le1\). By sign and permutation
we can also impose \(p=x_1+x_2\ge-x_7-x_8\); call this sorted half
polytope \(P\). In gaps \(g_i=x_i-x_{i+1}\) its nontrivial inequalities,
besides \(g_i\ge0\), are
\[
 2\sum_i g_i\le2,\qquad (1,2,3,4,3,2,1)\cdot g\le2,\qquad
 (-1,-2,-1,0,1,2,1)\cdot g\le0.
\]
Let \(h=(3/4,3/4,-1/4,\ldots,-1/4)\). Its coset modulo \(L\)
is invariant under signs and permutations: differences of permutations
are integral and zero-sum, and \(2h\in L\). It has squared torus distance
exactly \(3/2\), as the closest zero-sum vector with coordinates all
\(-1/4\) or all \(1/4\) modulo one has two coordinates of magnitude
\(3/4\) and six of magnitude \(1/4\).

The [small rational checker](verify.py) establishes four polytope facts.
Write \(q=\|x\|^2\) and
\(\gamma(u)=\max_{z\in Q}z\cdot u\), the gauge of \(R\).

1. The 21 vertices of \(P\) consist of \(h\), with \(q=3/2\), and 20
   vertices with \(q\le1\). In particular all squared torus distances
   are at most \(3/2\).
2. For each vertex \(y\ne h\), if the chord from \(h\) to \(y\) reaches
   \(\gamma(x-h)=1/2\), its first crossing has \(q\le1\).
3. On \(p\le3/4\) within \(P\), \(q\le7/8\).
4. On \(3/4\le p\le1\) within \(P\),
   \(q+2(1-p)^2\le1\).

For clarity, fact 2 is sufficient to show \(q\le1\) everywhere in
\(\{x\in P:\gamma(x-h)\ge1/2\}\). For any fixed \(z\in Q\), clip \(P\)
by \(z\cdot(x-h)\ge1/2\). All new vertices are on its chords; only chords
from \(h\) to a shallow vertex can matter. Their first crossings stay
no nearer to \(h\) than the first gauge crossing, so their norm bound
follows from convexity. Take the union over \(z\). Furthermore, any
\(x\ne h\) is of the form \(h+s(y-h)\) with \(0<s\le1\) and \(y\) a
convex combination of the shallow vertices. If \(q>1\), convexity
between \(h\) and the first gauge crossing (or between \(h\) and \(y\)
if there is no crossing) gives
\[
 \gamma(x-h)\le 3/2-q. \tag{1}
\]
The same holds without sorting, modulo the invariant coset \(h+L\).

Every point of \(R\) is at squared distance at most one from \(L\).
Indeed if \(x\in R\) violates \(x\cdot r\le1\), subtract \(r\). The
new vector remains in \(R\): for \(z\in Q\) with \(r\cdot z\ge0\)
the polar bound only improves; otherwise use
\(z-(r\cdot z)r\in Q\), the root reflection of \(z\), and \(x\cdot r>1\).
Root reflections preserve \(L\) and the roots because \(L\) is integral
and the roots have squared norm 2. Each subtraction strictly reduces
the norm; discreteness of the affine lattice in a bounded ball
forces termination at \(R\cap Q\), where polarity gives \(\|x\|^2\le1\).
Consequently any three torus phases have a pair of squared distance
at most one: if two displacements from one phase have squared distance
above one, both lie in \(h+R/2\bmod L\) by (1), so their difference
is represented in \(R\). The triple bound is sharp using \(0\) and two
half-integer shifts supported respectively on coordinates
\(\{1,2,5,8\}\) and \(\{1,3,4,8\}\), with two plus and two minus signs
on each support. Each of the three differences has four half-integer
coordinates, which remain four after either integral or half-integral
lattice shifts, and easily attains norm squared one.

## 2. The paired clearance estimate

For a torus offset \(x\) represented in \(P\), put
\(\sigma=\|x-h\|^2=q+3/2-2p\), an **upper** bound on the torus distance
to the coset \(h+L\). If \(q>1\), then \(\sigma<1/2\): in the same
decomposition \(x=h+s(y-h)\) with \(\|y\|^2\le1\), the exact squared
norm identity gives
\[
 q\le 3/2-\tfrac12s-\frac{1-s}{s}\sigma.
\]
If \(\sigma\ge1/2\), the right side is at most one by \(s+1/s\ge2\).
If \(7/8<q\le1\), facts 3 and 4 give
\(p\ge 1-\sqrt{(1-q)/2}\) (automatic when \(p\ge1\)).
Thus setting \(u=\sqrt{(1-q)/2}\in[0,1/4]\) gives
\[
q=1-2u^2,\qquad \sigma\le 1/2-2u^2+2u. \tag{2}
\]
For \(q>1\) the clearance bound in (2) remains true with \(u=0\).

## 3. The planar lattice cannot beat the incumbent

Take a Gauss-reduced normal basis \(A,B\), flipping a sign if necessary.
Write \(a=\|A\|^2\), \(c=\|B\|^2\), \(t=A\cdot B\), so
\(0\le2t\le a\le c\). The square of the normal covolume is \(ac-t^2\).
All normal steps have squared length at least \(1/2\) since the torus
clearance is at most \(3/2\). We compare the three fibers based at
\(0,A,B\), so no group or periodic rule for their horizontal offsets
is assumed.

**Both basis steps short, \(c<1\).** Their phases relative to the
fiber at zero are above one and therefore near the **same** \(h\) by (1).
Their difference is represented in \(X R\) where \(X=a+c-1\).
Every element of \(R\) has squared norm at most 2, so
\(a+c-2t\ge2-2X^2\). In particular \(X\ge1/2\) and
\(t\le(X-1/2)(X+1)\). Since \(ac\ge X\), if \(X\ge3/4\) then
\(ac-t^2\ge X-1/4\ge1/2\); otherwise
\(t^2\le(X-1/2)^2(X+1)^2\le X-1/2\), giving the same bound.

**Second step long, \(c\ge1\).** A candidate with \(ac-t^2<1/2\)
has \(a<2-\sqrt2\), since \(a\ge1\) already gives \(ac-t^2\ge3/4\),
and otherwise \(ac-t^2\ge a-a^2/4\). Write \(a=1/2+g\), \(0\le g<9/100\).
The phase of \(A\) has the form \(h+e\bmod L\), \(\gamma(e)\le g\).
For the phase \(x\) of \(B\), the torus distance squared to \(h+e\)
is at most \(\sigma+2g+2g^2\): use a nearest representative of \(x-h\)
in \(Q\), its dot-product bound with \(e\), and \(\|e\|^2\le2g^2\).
The sign and permutation selecting \(x\in P\) preserve the coset and
all necessary bounds.

If \(q\le7/8\), then \(c\ge9/8\) and
\(ac-t^2\ge9a/8-a^2/4\ge1/2\). Otherwise use (2); together with \(c\ge1\)
this gives \(c\ge1+2u^2\) and, from the \(B-A\) step,
\[
c\ge1+2u^2+2(t-u)-(3g+2g^2).
\]
Write \(h_0=(3g+2g^2)/2\). The lower bound
\(a(1+2u^2)+2a\max(0,t-u-h_0)-t^2\) decreases to the breakpoint
and then increases on \(0\le t\le a/2\). Its minimum is at least
\(a(1+2u^2)-(u+h_0)^2\). This last expression decreases in
\(u\in[0,1/4]\), because \(4gu-2h_0\le0\). Evaluating there gives
\[
ac-t^2\ \ge\ \frac12+\frac38g-
       \bigl(\tfrac12+(\tfrac32+g)^2\bigr)g^2 \ \ge\ \frac12.
\]
So the planar covolume is at least \(1/\sqrt2\). The full covolume
per point is at least \(1\), since the horizontal covolume is \(\sqrt2\).
Radius \(1/\sqrt2\) now yields the claimed center density.

Equality is possible: take perpendicular normal steps with lengths
\(1/\sqrt2,1\), assign additive phases \(h,x\), respectively, where
\(x=(1,1,-1,-1,0,0,0,0)/2\).
The three relevant squared torus clearances are \(3/2,1,1/2\) for
\(h,x,x\pm h\). All farther normal displacements already have squared
length at least 2. This recovers a tie, not a new density.

References: [Conway–Sloane, *Laminated lattices* (1982)](https://annals.math.princeton.edu/1982/116-3/p07),
[catalogue with exact \(\Lambda_9\) comparison](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html),
[project density table](https://github.com/alexlegeartis/SpherePackings/blob/main/LOWER_BOUNDS_1_96.md),
and [Conway–Sloane's classical E8-fiber and fluid-diamond ties](https://neilsloane.com/doc/BEST/Me194.pdf).

[Extension to squared minimum in [29/25, 2]](RADIUS_EXTENSION.md).
