# Local optimality of the three-translate 55D packing

The three-translate packing of center density **92,197.2234049774…** is a
strict local density maximum among 55D packings with three translates per
period. This allows simultaneous changes to every period vector, all
motif coordinates and the sphere radius. The statement is up to Euclidean
isometries and overall scale. [Construction and minimum](README.md).

The proof enumerates the **2,499,741** closest points in the 54D coset
$Q+L$. Each short gap has that many contacts per center with the next
layer. A long gap has **584**. The long-gap contacts span the directions
needed for strictness. Their exact moments and a positive stress on 8,748
parent roots certify the local result.

## Enumerating the closest points

The coefficient field is $\mathbb F_{3^6}$ and the curve is
$y^2=x^3+x+t^{28}$. The source height formula for a reduced fraction
$x=a/b$ is
\[
\widehat h(x,y)=\max(\deg a,28/3+\deg b).
\]
At height $D=28/3$, $x$ is a polynomial of degree at most 9. Hence
$y$ has degree 14 with leading coefficient $1$ or $-1$. Work with
leading coefficient $1$ by taking the negative point if necessary.
Translation of $t$ by a unique field element normalizes the coefficient
of $t^{13}$ in $y$ to zero. Each translation lifts to the elliptic curve
by an additive correction to $x$. [Published height and rank inputs](https://doi.org/10.1007/s40993-022-00321-0).

Let $A=y^2-t^{28}$. In degrees above 9, its coefficients at exponents
indivisible by 3 vanish since $A=x^3+x$. Comparing those coefficients
from the top gives
\[
y=t^{14}+u t^{10}+v t^7+u^2t^6+w t^4-uv t^3-u^3t^2+zt+v^2-uw.
\]
The equation $A_{11}=0$ gives $uz+vw=0$. Using this equation in
$A_{21}=A_7^3$ gives $v^9=v$. For $u\ne0$, the equation
$A_{24}=A_8^3$ is the quadratic
\[
(w^3)^2+(v^6/u^3)w^3+u^{15}+u=0,\qquad z=-vw/u.
\]
For $u=0$, the same equations give $w=z=0$. In detail, $vw=0$ gives
$w=0$ when $v\ne0$, then $A_{24}=A_8^3$ gives $z=0$.
When $v=0$, that equation gives $w=0$, then $A_{15}=A_5^3$ gives
$z=0$.

The verifier checks every root of these quadratics for all $u$ and
$v\in\mathbb F_9$. It reconstructs the nonconstant coefficients of
$x$ from $A_j=x_j+x_{j/3}^3$, omitting the second term for $3\nmid j$.
It then checks every remaining coefficient and solves $x_0^3+x_0=A_0$.
There are **1,143** normalized values of $y$, each with three choices of $x_0$.
One choice for each $y$ is bound to $Q+L$ using its exact pairings
with a full-rank set of generators. The generator isometries for
$x\mapsto x+\sqrt{-1}$ and the six field-basis translations preserve
$Q+L$. Applying them gives $1,143\cdot3\cdot729$ distinct points in
that coset. Negation gives the complete shell in $-Q+L$.
[Parent and coset reconstruction](../cross_coset55_20261116/PROOF.md).

A closest point $m$ contributes to a long gap precisely when
$\langle m,Q\rangle=-D/2$. Before undoing the normalization, the
verifier computes this inner product by applying the inverse translation
to $Q$ instead of to $m$. This tests the full shell with a small
matrix product. The contact set $C$ has 584 points. The involution
$m\mapsto-m-Q$ preserves $C$, so $z_m=m+Q/2$ has average zero and
lies in $Q^\perp$. The checked rank of these $z_m$ is 53.

## Strict local density maximum

Fix squared minimum at 10 by rescaling. Near the given packing one can
use the first 54 period generators to define the horizontal space.
After an isometry, a nearby rank-55 period with a three-point motif
is described by a perturbed horizontal Gram $G$ and three successive
gaps. Each gap has a horizontal offset and a vertical length. Let $B$
be the first-order change in $G$ and $s_i$ the change in the horizontal
offset of gap $i$. Let $u_i$ be the change in its vertical length measured
in units of $h=\sqrt{2/3}$.

Write $a=1/5$, $g=(1,1,9/5)$ and $H=19/5$. The first variation
of the logarithm of the period covolume is
\[
\frac12\operatorname{tr}(G^{-1}B)+\frac{u_1+u_2+u_3}{H}.
\]
For a horizontal contact $v$ across gap $i$, the first variation of
its squared length is $v Bv^T+2\langle v,s_i\rangle+2h^2g_i u_i$.
For the period contact use $v=aQ$, the sum of gap changes, and $H$.
For a parent contact use $v Bv^T$. Each such variation is nonnegative
along a feasible first-order displacement.

Use the following weights on these variations. A group's weight is
distributed uniformly over the indicated set. Put $p=3/76$.

| Contact group | Total weight |
|---|---:|
| Period contact | $p$ |
| For each short gap, the 2,187-point isotropic subset of $Q+L$ | $4p/5$ |
| For each short gap, the full set $C$ at offset $Q$ | $2p/5$ |
| Long gap, the full set $C+aQ$ | $2p/3$ |

The isotropic set has mean zero and second moment $(D/54)G^{-1}$
in basis coordinates. The other two sets have means $-Q/2$ and
$-3Q/10$. Thus the weights cancel the $s_i$ contributions and give
coefficient $1/H$ for each $u_i$. Let $Z$ be the mean of $z_m^Tz_m$
over $m\in C$. The remaining horizontal tensor needed to complete the
log-covolume variation is
\[
S=\left(\frac12-\frac{8pD}{5\cdot54}\right)G^{-1}
  -p\left(\frac{22}{15}Z+\frac3{10}Q^TQ\right).
\]
The verifier expresses $S$ as a **strictly positive** combination of
the 8,748 parent root tensors. It solves one rational linear system
from the sparse stress seed and checks the tensor identity and all
signs exactly. [Stress seed](stress_seed.json).

If a first-order displacement is feasible and has nonincreasing
covolume, all positively weighted contact variations must vanish.
The parent root tensors span the symmetric matrices, so $B=0$.
For each short gap the mean-zero contact subset forces $u_i=0$.
Its full span then forces $s_i=0$. For the long gap, subtracting the
variations of the paired contacts $m+aQ$ and $-m-Q+aQ$ makes $s_3$
orthogonal to all $z_m$. Thus $s_3$ is parallel to $Q$.
The long-gap and period averages give two independent equations in
$\langle Q,s_3\rangle$ and $u_3$, forcing both to vanish. There is no
nonzero feasible first-order displacement with nonincreasing covolume.
If better or equal noncongruent packings approached this one, a limit
of their normalized parameter differences would give such a displacement.
This proves strict local optimality.

The exact checks use NumPy and python-flint.

~~~sh
python research/uneven_layers55_20261122/verify.py
python research/uneven_layers55_20261122/verify_local_optimum.py
~~~
