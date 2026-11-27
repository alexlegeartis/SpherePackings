# The extremal even 48D stratum has no two-neighbor edges

Two distinct even unimodular dimension-48 lattices of minimum 6
cannot have intersection of index two in both. This exact
arithmetic obstruction helps explain why enlarging a canonical
normal subgroup does not give some seemingly large 47D gains.

Suppose \(L_0,L_1\) were such a pair, and set \(K=L_0\cap L_1\).
Then \(K\) is even, has minimum at least 6, and
\(K^*/K\) is a nondegenerate group of order four. The two different
subgroups \(L_i/K\) each have order two and are isotropic for the
quadratic form \(\|x\|^2/2\bmod1\). Consequently the discriminant
group is the Klein four group; the two nonzero even classes must
pair by \(1/2\bmod1\), by nondegeneracy. Their sum therefore has
quadratic value \(1/2\bmod1\). Adjoining this third nonzero class
gives an **odd** integral unimodular lattice \(U\), whose even
sublattice is precisely \(K\). Its shadow
\(K^*\setminus U\) is the union of the two even nonzero classes and
has minimum at least 6.

Use \(q=e^{\pi i\tau}\), and write the standard odd-unimodular theta
expansion as
\[
 \Theta_U=\sum_{j=0}^6 a_j\theta_3^{48-8j}\Delta_8^j,\qquad
 \Theta_{\rm shadow}=\sum_{j=0}^6\frac{(-1)^j a_j}{16^j}
         \theta_2^{48-8j}\theta_4(2\tau)^{8j},
 \quad \Delta_8=\theta_2^4\theta_4^4/16.
\]
The leading shadow exponents for \(j=6,5,4\) are respectively 0, 2
and 4, so successively \(a_6=a_5=a_4=0\).
The coefficient \(a_0=1\). There is no norm-one vector in \(U\)
(its double would have norm four in \(K\)), and the even sublattice
excludes norms two and four. Hence the \(q,q^2,q^4\) coefficients
of \(\Theta_U\) vanish. The elementary expansion gives
\(a_1=-96,\ a_2=2400,\ a_3=-12700\). It follows that \(\Theta_U\)
has *exactly 100 norm-three vectors*. See the [integer replay](verify_shadow.py).

But any two distinct antipodal pairs of norm-three vectors are
orthogonal: their sums and differences are nonzero vectors of \(K\),
so \(6\pm2(x\cdot y)\ge6\). Hence there are at most \(48\) antipodal
pairs, or \(96\) such vectors. This is a contradiction.

## Consequence for the 47D order-two architecture

In the [normal-slice construction](README.md), take a primitive
normal \(v\) of squared norm \(4t\), where \(t\in\{2,4,6,8,10\}\).
The two-sheet subgroup consists of the zero and opposite classes
modulo \(4t\). Its even unimodular 48D neighbor
\(\Lambda'=K+\mathbb Z(v/2)\) is distinct from the original by
primitivity, so the result above forces minimum at most 4 in
\(\Lambda'\). Write \(r=v/2\), of squared norm \(t\).
The same weight-26 harmonic test from the main note shows \(r^\perp\)
in \(\Lambda'\) has minimum at most 4 for all the stated \(t\).
Indeed if there are roots but none perpendicular, the \(q\)
harmonic sum is positive (each contribution at least \(1-t/24\))
and the \(q^2\) coefficient is negative, impossible unless some
norm-four vector is perpendicular (otherwise each contribution
is at least \(1-t/12>0\)). If there are no roots, the same harmonic
series is zero and the lattice has a nonempty norm-four shell.
Thus the opposite phase has squared clearance at most 4; the
zero phase alone has inherited minimum at least 6.
For \(t=2,4\) the more general \(b\)-neighbor of the main note
does not need the restriction \(b=2\).

This is an obstruction for **subgroup** projections, not a claim
about arbitrary pairs of displacements for larger normal lengths,
nonintegral regroupings, or absolute 47D optimality.

The theta and shadow expansions are the standard formulas of
[Rains–Sloane, *The Shadow Theory of Modular and Unimodular Lattices*](https://arxiv.org/abs/math/0207294);
for [48D neighbors from ternary codes, see Nebe, §3](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/DIM48.pdf).
