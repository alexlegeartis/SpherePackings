# Seven complete 88D golden trace curves

For **all seven** cyclotomic ideal lattices over \(\mathbb Q(\zeta_{115})\)
in the 88D section of the [Nebe–Sloane catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/unimodular.html),
we exhibit literal squared-norm-8 vectors of both pairing signs \(\pm8\).
The expected gain on this family therefore fails: **the full golden
trace curve has exact Hermite optimum 8**, attained only by the
unwarped parent up to scaling and golden-unit equivalence. In
particular the 5-modular endpoints have exact squared minimum 16,
not 20. This is a scoped result, not any bound on arbitrary 88D
lattices or packings.

| Parent | Source | Endpoint witness \(q,s\) | \(\min(xMx^\top)\) |
| --- | --- | --- | ---: |
| ideal115L0 | [catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/dim88a.html) | \(8,-8\) | 16 |
| ideal115L | [catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/dim88b.html) | \(8,-8\) | 16 |
| X139 | [catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/dim88c.html) | \(8,-8\) | 16 |
| X599 | [catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/dim88d.html) | \(8,-8\) | 16 |
| X4049 | [catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/dim88e.html) | \(8,-8\) | 16 |
| X4969 | [catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/dim88f.html) | \(8,-8\) | 16 |
| X1151s24 | [catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/dim88g.html) | \(8,-8\) | 16 |

## Exact verification

Let \(G\) be the listed Gram matrix, and \(U\) the first automorphism,
expressed in the row convention \(UGU^\top=G\). Set
\[
 Q=U^{46},\quad S=I+2(Q+Q^4),\quad
 q(x)=xGx^\top,\quad s(x)=xSGx^\top,\quad M=(5G+SG)/2.
\]
The standard-library verifier checks positive definiteness and
determinant one of every \(G\), the displayed automorphism equation,
\(Q^5=I\), \(\sum_{i=0}^4Q^i=0\), \(S^2=5I\),
\(SG=(SG)^\top\), and the exact determinant \(5^{44}\) of \(M\).
The golden unit \(T=(I+S)/2\) is integral and invertible over
\(\mathbb Z\), with \(T^{-1}=T-I\).
The file witnesses.json gives one vector of norm 8 and pairing
\(-8\) in each published basis. Applying \(T\) sends its \((q,s)\) to
\[
  ((3q+s)/2,(5q+3s)/2)=(8,8),
\]
and that second literal vector is also verified. Discovery used
floating-point block reduction, not relied upon in the proof.

Run python -B verify.py in this directory. The file inputs.bin.gz
is a copy of the mathematical matrix inputs from the linked catalogues
and [Gauvin's source archive](https://github.com/Katabolanga/rang88-reseaux-ideaux).
The extraction retains only the lower triangle of \(G\) and the full
\(U\), in row-major order as signed little-endian int16, for each
parent in the order of the table. Its compressed SHA-256 is
8741db7e6453bfeb445b3a5722359745ab1145be3e44985f0a51b42d2f4846c1.
There are no numerical solver inputs or heuristic inequalities in
the exact replay.

## Exact curve and scope

Every positive metric in the golden trace family, modulo scale and
integral pullback by powers of \(T\), is \(G+tSG\) with \(|t|\le1/5\).
In fact the ratio of its two embedding weights is
\((1+\sqrt5 t)/(1-\sqrt5 t)\), and pullback by \(T\) multiplies that
ratio by \(((1+\sqrt5)/2)^4\). The two endpoints differ by this factor.
Each \(S\) eigenspace has rank 44 (irrational conjugate eigenvalues of
an integral matrix), so the normalized Hermite denominator is
\(\sqrt{1-5t^2}\).

For \(t\ge0\), the literal \((8,-8)\) vector gives a squared minimum
at most \(8(1-t)\); for \(t\le0\), use the literal \((8,8)\) vector.
For the matching lower bound, use the original minimum 8 proven
or reported in each catalogue's cited primary construction, not
reproved by this verifier. Every nonzero \(x\) then satisfies
\(q(x)\ge8\) and \(q(xT^{\pm1})=(3q\pm s)/2\ge8\).
Consequently
\[
 q+ts\ge (1-3|t|)q+16|t|\ge8(1-|t|),
\]
which gives the exact stated minimum curve. On the fundamental interval,
\[
 \gamma(t)=\frac{8(1-|t|)}{\sqrt{1-5t^2}}\le8
\]
with strict inequality for nonzero \(t\), since
\(6|t|^2-2|t|<0\) there. In particular the upper bound of 8 for
the full trace family does **not** depend on a separate proof
of the published original minima; those are used only for sharp
attainment and the reverse inequality.

The resulting packing fraction is less than \(4.959\cdot10^{-20}\)
using \(\pi<355/113\), under 26.8% of the pinned
[88D lattice benchmark](https://github.com/alexlegeartis/SpherePackings/blob/6696a479a8a55e440ff8e65f10702702fe7d599c/README.md).
This proves no upper bound beyond the named trace families; in particular
the other golden rank-88 parents or non-trace deformations are not covered.
