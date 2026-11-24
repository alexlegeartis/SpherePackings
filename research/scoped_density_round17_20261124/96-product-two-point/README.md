# All-degree two-point barrier for the balanced norm-12 golden shell

In dimension 96 let \(L\) be even unimodular with minimum at least 8,
and let \(U\) be an integral self-adjoint operator satisfying \(U^2=U+1\).
Write \(S=2U-I\), \(q(x)=x\cdot x\), \(s(x)=x\cdot Sx\), and assume
\(q(x)=8\) implies \(s(x)=0\). Let \(A\) count norm-8 vectors.
The [previous trace calculation](../../scoped_density_obstructions_round15_20261123/96-golden-trace/README.md)
gives \(0\le A\le1,662,897,600\) and
\(N=1,218,201,435,360,000+3,212A\) vectors in the balanced norm-12 shell.
The [split-shell design result](../../scoped_density_obstructions_round16_20261124/96-golden-shells/HILBERT.md)
supplies both projected covariance identities and the mixed covariance
identity on this balanced shell.

**Result.** For *every* such \(A\), the single-shell, product-sphere
two-point relaxation described below has a formal rational distance
distribution with exactly that \(N\). It satisfies all degrees of
the product Gegenbauer inequalities, the known second moments,
the pair-span norm and shell tests through norm 16, and the
corresponding shell-size bounds for each separate pairing label. It
even obeys a maximum-dot degree budget of **753,000 per sign**, more
than 14 times below the proved budget of 11,032,662.
Increasing the polynomial degree, or tightening just that local
budget to 753,000, therefore cannot exclude this branch with those
two-point inputs. This is not a point configuration, a lattice, or
a new upper or lower bound on density. It does not rule out stronger
relaxations coupling different shells or using three-point data.

## The finite relaxation and witnesses

For an ordered pair of nonidentical, nonantipodal balanced norm-12
vectors \(x,y\), use \(p=x\cdot y\) and \(b=x\cdot Sy\). Each eigenspace
has dimension 48 and the projections have norm \(\sqrt{6}\); the two
normalized dots are \(t_\pm=(p\pm b/\sqrt5)/12\). The averaged distance
distribution divides ordered pair counts by \(N\). Antipodality and
averaging under the formal operation \(S\mapsto -S\) permit grouping
the values with opposite signs of either \(p\) or \(b\). This averaging
does not presume an isometry swapping the two eigenspaces. Count
entries \(w_{p,b}\) in the certificate are the total weights on these
orbits, for \(p,b\ge0\). The diagonal and antipodal pair each have
weight 1. The design equalities are
\[
  2+\sum w=N,\qquad
  288+\sum w p^2 =3N/2,\qquad
  \sum w b^2=15N/2.
\]
The mixed \(pb\) moment vanishes after sign averaging. Equivalently
the pure factor degree-2 and the mixed \((1,1)\) Gegenbauer moments
are zero. All odd total degrees vanish by antipodality.

For each supported \(p,b\), the checker builds the integral Gram of
\(x,Ux,y,Uy\). It is positive definite (positive elimination
pivots), and inverse-Gram coordinate bounds enumerate all vectors of
squared norm at most 16 in that integral span; in particular none
have \(q<8\). For any such vector with nonzero coefficient pair on
\(y,Uy\), fix \(x\) and the signed label. The map from \(y\) into
its resulting signed shell is one-to-one, since every nonzero
\(\beta\in\mathbb Z[U]\) acts invertibly. Thus an individual label has
at most the corresponding shell population; the checker tests all
these bounds using the exact Hilbert theta counts for \(q\le16\).
Signs conjugate the reached shell or leave it unchanged, and those
counts are equal for both signs. In particular no positive-weight
label forces a vector in a missing shell, including norm 8 when
\(A=0\) or signed norm 10 with \(|s|=6\) at the other endpoint.
The entries are averaged counts and need not be integers.

These constraints are affine in \((A,w)\), including the per-angle
injectivity bounds and the max-dot budget. The two rational lists
in the counts file certify the endpoints; affine interpolation covers
the entire interval without any assumption of intermediate lattice
existence. The two average counts with \(p=8,b=0\) per sign are 0
and about 752,865, respectively. The support test is independent
of \(A\), so interpolation does not introduce new labels.

## All degrees at once

Let \(P_k\) be the normalized Gegenbauer polynomial for \(S^{47}\).
For even total \(i+j\) the harmonic requirement is
\[
   2+\sum_{p,b} w_{p,b}
   \frac{P_i(t_+)P_j(t_-)+P_i(t_-)P_j(t_+)}2 \ \ge\ 0 .
\]
This is the sum-of-squares inequality from the spherical addition
formula and the product of the two zonal kernels. The checker
evaluates both endpoint lists for \(i,j<100\) by outward-rounded
integer intervals. Apart from the exact zero design modes, every
finite tested nonconstant moment exceeds 0.045.

For completeness the tail does not require a high-degree solver.
The [Laplace integral for the normalized Gegenbauer polynomial](https://dlmf.nist.gov/18.10.E4)
writes \(P_k(t)\) as the expectation of
\((t+\mathrm{i}u\sqrt{1-t^2})^k\), where \(u\) has density
proportional to \((1-u^2)^{22}\) on \([-1,1]\). The magnitude of the
base is at most 1. For every \(k\ge100\),
\[
  |P_k(t)|\le
  E[(t^2+(1-t^2)u^2)^{50}].
\]
The even moments of \(u\) are obtained successively by multiplying
by \((2j-1)/(2j+45)\). Taking the other factor at most 1 in absolute
value, the checker sums the expectation bound at both embeddings
with the exact weights. At either endpoint the whole off-diagonal
tail is less than 0.087, below the diagonal-plus-antipodal
contribution 2. All Gegenbauer modes are thus verified by one
finite calculation; no numerical solver result is used as proof.

Run **python -B replay.py** and **python -B theta.py** (standard
library only). The latter extends the scalar Hilbert theta
calculation through squared norm 16 and independently matches each
signed coefficient used by the verifier. The full parallel
weight-24 space is the weight-24 monomials of the signed symmetric
ring \(\mathbb C[G_2,G_5,G_6]\), together with the weight-9
monomials times the weight-15 generator \(G_{15}\). This uses the
[previous Aoki–Takemori ring and generator calculation](../../scoped_density_obstructions_round16_20261124/96-golden-shells/HILBERT.md).
Imposing scalar constant 1, zero coefficients below squared norm 8,
and only the norm-8 balanced coefficient \(A\) gives full rank for
the space. The script solves separately for its constant and
\(A\)-linear parts. The covariance input is the already proved
design identity, not an assumption of realizability.
