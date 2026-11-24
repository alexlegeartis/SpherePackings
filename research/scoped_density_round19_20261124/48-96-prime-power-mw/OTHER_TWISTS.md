# The other tame CM twists do not rescue the prime-square tower

The [same integral argument](README.md) applies to the quadratic
and cubic twists in \(p=5,q=25\) and the quadratic twist in \(p=7,q=49\).
These exhaust the nontrivial standard exponent quotients of the tame
CM groups on the respective supersingular elliptic curves. This is
only the base polynomial \(t^q-t\), not arbitrary superelliptic bases
or sublattices. The unique supersingular elliptic \(j\)-invariants in
characteristics 5 and 7 here are \(0\) and \(1728\), respectively:
the Hasse invariant is the coefficient of \(X^{p-1}\) in
\((X^3+AX+B)^{(p-1)/2}\), respectively \(2A\) and \(3B\).

| Rank | Equation, \(f=t^q-t\) | \(\tfrac12 v_p(\det\mathrm{MW})\) | Squared min upper | Center-density ratio cap |
|---:|---|---:|---:|---:|
| 48 | \(y^2=x^3+f^3\) | 24 | 10 | \((2/3)^{24}\) |
| 48 | \(y^2=x^3+f^2\) | 24 | \(40/3\) | \((16/27)^{12}\) |
| 96 | \(y^2=x^3+f^2x\) | 48 | 14 | \((2/9)^{24}\) |

The references are \(P_{48}\) and Bos–Conway–Sloane 96D respectively.
In the crystalline calculation keep the two distinct elliptic Hodge
lines indexed by \(b=1,-1\), even if both restrict to the same
character modulo 2. Their Frobenius exchanges them (in both
characteristics, acting on the full order 6 or 4 group).
The active pairs have \(a\equiv b\bmod k\) and
\((a,b)\mapsto(pa,-b)\). The exponent on each pair is
\(\epsilon(a)+\mathbf1_{b=-1}-1\).
Run verify.py for the exact orbit sums and density arithmetic.

For the quadratic cases the finite and infinite fibres are
all \(I_0^*\), respectively 26 and 50, with root determinant
4 each. Torsion consists of the four fixed constants \(E_0[2]\).
Thus the covolume lower bounds from Shioda's formula are
\(5^{24}/2^{24}\) and \(7^{48}/2^{48}\).
For the cubic rank-48 case there are 25 type IV and one type IV*
fibre, each of root determinant 3, and the fixed constants
have order 3. This gives covolume at least \(5^{24}/3^{12}\).

## Three low-degree elliptic quotients before the prime-square pullback

Use \(s=t^p+t\), which gives \(s^p-s=t^{p^2}-t\)
by a separable degree-\(p\) base change.
After passing from the degree-\(p\) twist to
\(C:v^k=s^p-s\), a degree-\(d\) equivariant map to
\(E_0\) gives a section of squared height \(2d/k\)
before that degree-\(p\) pullback. Indeed on the constant
curve its height is twice the intersection with the zero
section (degree \(d\)); height scales by the base extension
degree. The unique infinity of \(C\) is fixed by the tame
action; send its image to the origin of the elliptic quotient.

* **\(p=5,k=2\)**. On the genus-two curve take \(s\mapsto
  a/s,\ v\mapsto b v/s^3\) with \(a=2,\ b^2=-a\).
  For \(a\in\mathbb F_5^*\) the new value of \(s^5-s\) is
  \(-a(s^5-s)/s^6\), and \(a^2=-1\) makes the lift an involution.
  The two fixed \(s\)'s are opposite roots of \(s^2=a\);
  neither is a branch value. The \(v\) multiplier is \(+1\)
  at one and \(-1\) at the other, so there are two fixed
  points. Riemann–Hurwitz gives an elliptic quotient of
  degree 2. It is supersingular (the source Jacobian is
  supersingular by the same orbit test at \(q=p\)).
  The hyperelliptic involution induces negation on the
  quotient, since it acts by negation on the Jacobian.
* **\(p=5,k=3\)**. On the genus-four curve take the commuting
  involutions \((s,v)\mapsto(-s,-v)\) and \((s,v)\mapsto
  (1/s,-v/s^2)\). The first fixes \(0,\infty\), the second
  fixes the two branch points over \(s=\pm1\), and
  their product fixes the two branch points over \(s^2=-1\).
  They have disjoint fixed sets. Riemann–Hurwitz gives
  genus one for the degree-four quotient. The deck
  \(\mu_3\) commutes with these maps and acts on its
  differential by \(\zeta^{-1}\), seen on the invariant
  differential \(ds/v\). This is the needed orientation
  of the action on the supersingular \(j=0\) curve.
* **\(p=7,k=2\)**. On the genus-three curve take
  \((s,v)\mapsto(-1/s,v/s^4)\). It is an involution:
  \((s^7-s)\) transforms to \((s^7-s)/s^8\).
  The two fixed \(s\)'s satisfy \(s^2=-1\), are not branch
  values, and each fixes two \(v\)'s. The degree-two
  quotient has genus one and is supersingular.
  The deck involution acts by negation as in the first case.

Multiplying the resulting heights by \(p\) gives the stated
short vectors. The regulator cost remains the obstacle for all
three twists even with these local component corrections.
