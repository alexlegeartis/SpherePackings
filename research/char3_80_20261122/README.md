# An 80D Mordell–Weil lattice and a small-index obstruction

The curve $y^2=x^3+x+t^{41}$ over $\mathbb F_{3^{16}}(t)$ gives an
explicit rank-80 height lattice of Gram determinant $3^{48}$ and
squared minimum 14, which gives 1.68% of MW80's density. At minimum 16
it would beat MW80 by a factor of about 3.51. Even an index-three
sublattice of minimum 16 would beat MW80 by 17.04%. Every index-two and
index-three sublattice still has minimum 14 by the obstruction below.
[MW80 comparison and source reconstruction](../creative_sections_20261116/PROOF.md#78d-and-79d-antipodes-from-the-reconstructed-mw80).

## Exact lattice and minimum

We found the initial sections by pulling the curve back along $t=u^2$
and tracing polynomial points under $u\mapsto-u$. The verifier starts
from the resulting sections. It checks all their curve equations and
height pairings, the exact span and its determinant $3^{80}$.
The next 16 generators come from signed sums divided by the
degree-three endomorphism $\phi=\rho-1$, where
$\rho(x,y)=(x+\tau,y)$ and $\tau^2=-1$. In affine coordinates,
\[
\phi(x,y)=(-x^3-t^{41}-\tau,\ \tau y^3).
\]
Taking cube roots verifies each division in the function field. Exact
Hermite normal form of all the generators gives the final determinant
$3^{48}$. [Certificate](certificate.json).

The norm of a point is $2\deg(y)/3$, where $\deg(y)$ is the degree
of the rational function on the line. Indeed, the Artin–Schreier cover
$v^3+v=t^{41}$ has degree three, and $(x+v,y)$ is a map to the
constant curve $y^2=X^3+X$. The degree of this map is $\deg(y)$,
since the coordinate $y$ on the constant curve has degree three.
The $X$ coordinate has degree two. The degree of multiplication
by $m$ on the constant curve is $m^2$, and the fixed change
$X=x+v$ alters the coordinate height by a bounded amount. This gives
the stated canonical norm after dividing by the cover degree.
The resulting form is positive on nonzero points. A constant image
would be fixed by $X\mapsto X+\tau$, whose only fixed point is the
identity. Polarization gives the pairings used by the verifier.

For any nonzero point, $x$ must have a pole at infinity of order at
least 14. Otherwise $t^{41}$ would strictly dominate the right-hand
side and the pole order of the square $y^2$ would be odd. With a pole
of order at least 14, $x^3$ dominates and $y$ has pole order at least
21, giving squared norm at least 14. The verifier explicitly binds two
polynomial points of degree $\deg(y)=21$ to the final lattice.
The squared minimum is therefore 14.

## Short vectors in every index-two and index-three sublattice

For a norm-14 vector $v$, all of $v,\rho v,\rho^2v$ have norm 14 and
sum to zero. A homomorphism to $\mathbb F_2$ must vanish on at least one
of them. Thus every index-two sublattice contains a norm-14 vector.

For a homomorphism $\ell$ to $\mathbb F_3$ to avoid norm-14 vectors in
its kernel, its values on each such triple must all be the same. In
particular, $\ell(\phi v)=0$ on the norm-14 vectors. Since
$\phi^2=-3\rho$, the quotient by $\phi L$ is a vector space over
$\mathbb F_3$. The verifier computes its dimension as 40.

The two short witnesses generate **1,312** distinct projective rays in
this quotient using $t\mapsto\zeta t$ with $\zeta^{41}=1$, coefficient
Frobenius, and $(x(t),y(t))\mapsto(-x(-t),\tau y(-t))$. The verifier
checks the exact integral isometries and that the rays span the
quotient. Multiplication by $\phi$ identifies $L/\phi L$ with
$\phi L/3L$, so $\ell$ must factor through $L/\phi L$.

For each ray $e$, any allowed functional would satisfy $\ell(e)^2=1$.
The vector of these squares lies in the column span of the matrix of
all quadratic monomials $e_i e_j$. That matrix has rank 808 over
$\mathbb F_3$. Appending the all-ones column raises its rank to 809.
This rules out the functional and proves the index-three obstruction.

At squared minimum 16, index four or greater would already cost too
much covolume to exceed MW80. The calculation concerns this specific
lattice and these indices. Other curves, changes to the metric, and
other periodic motifs are outside its scope.

## Replay

With NumPy and python-flint, run

~~~sh
python research/char3_80_20261122/verify.py
~~~

The arithmetic uses polynomial operations over finite fields and exact
matrix operations. The density comparison uses rational squares.
