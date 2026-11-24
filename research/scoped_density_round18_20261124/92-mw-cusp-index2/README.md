# MW92: all index-\(47^2\) sublattices through the second cusp jet are excluded

Use the 92-dimensional Mordell–Weil lattice \(L\) of
\(Y^2=X^3+t^{47}-t\) and translation \(T:t\mapsto t+1\).
Write \(D=1-T\), \(K_j=D^jL\).
For **every** sublattice \(M\supset K_2\) of index \(47^2\) in \(L\),
the squared minimum is at most 18, and \(K_1\) attains 18.
Thus none beats the parent of squared minimum 16 and covolume \(47^7\):
the maximal center-density ratio in this family is
\((18/16)^{46}/47^2<1\).

In particular this excludes every index-\(47^2\) **translation-invariant**
sublattice without imposing the extra jet condition. Indeed completion
at 47 gives the action of \(\mathbb Z_{47}[\zeta_{47}]\) with \(D\) a
uniformizer and residue field \(\mathbb F_{47}\)
(\(\Phi_{47}(1-X)\) is Eisenstein).
A quotient of order \(47^2\) has module length 2, hence is annihilated
by \(D^2\). Completion at 47 does not change this finite quotient, so
\(K_2\subset M\).

The [previous exact reduction](../../scoped_density_round17_20261124/92-MW-cusp-jets/README.md)
left six isometry representatives of planes in the two-jet quotient.
Here are their literal short-vector squared norms:
\[
\begin{array}{c|rrrrrr}
(a,b,c)&(1,1,0)&(1,2,1)&(1,2,27)&(1,3,0)&(1,3,12)&(1,6,24)\\
\hline
\hbox{squared norm}&16&16&16&16&18&18.
\end{array}
\]
The six integer rows are in [witnesses.json](witnesses.json).
In the same checked quotient \(v\mapsto(u_1,u_2,j_1,j_2)\),
a row must have \(u_1\ne0\), \(u_2=b u_1\), and
\(u_1j_2-u_2j_1=c u_1^2\) modulo 47.
These are exact plane-membership tests, not a claim that the sampled
short vectors exhaust any shell.
The prior checker also handles all planes whose projection is not a
parallel one-dimensional line, so this is exhaustive for the stated index
and containment. It does not exclude arbitrary larger-index sublattices,
non-invariant sublattices outside this quotient, or different lattice
constructions.

Replay with OPENBLAS_NUM_THREADS=1 python verify.py (Python 3 and NumPy).
The previous checker is run again and supplies the exact verified Gram,
quotient and orbit data; its path can also be passed as the first argument.

## Partial next index: \(47^3\)

The same computation excludes improvement by **26,556 of 106,080**
translation-invariant sublattices of index \(47^3\).
Over the completed discrete valuation ring the two elementary-divisor
profiles are \((1,2)\) and \((0,3)\), in \(D\)-adic exponents. There
are 48 of profile \((1,2)\), one for each line strictly between \(K_2\)
and \(K_1\) of the required index. They contain a norm-20 vector of \(K_2\).

Those of profile \((0,3)\) are primitive lines modulo \(D^3\), a
projective line over the length-three local ring. Each of the 48
projective lines modulo \(D\) has \(47^2\) lifts. For any such
sublattice with leading residue line \(l\), its intersection with
\(K_2/K_3\) is exactly \(D^2l\) under the natural identification
induced by \(D^2\). Thus a norm-20 witness in that secondary residue
line rules out every lift. The original norm-20 witness and one new
literal norm-20 row, together with the three checked isometries
normalizing \(T\), give twelve different leading residue lines.
This excludes \(48+12\cdot47^2\) sublattices because
\(20^{46}<16^{46}47^3\). The other leading directions remain open
for this compact certificate. No assertion about all index-\(47^3\)
sublattices is made.
