# 47D: complete projection-phase ceilings for five further normal lengths

**Scope.** Let \(\Lambda\) be any extremal even unimodular lattice of
rank 48 (squared minimum 6). Let \(v\in\Lambda\) be primitive with
\(w=\|v\|^2\in\{12,14,16,18,20\}\). Over \(L=\Lambda\cap v^\perp\),
allow any union of translates by classes in **all** of \(L^*/L\).
At every possible radius, its center density is strictly below the
[Chen–Hu–Li–Wang–Wu 47D construction](https://arxiv.org/html/2505.02394v1#S3.SSx5),
\(\delta_*=47^{47/2}2^{-118}\simeq5925.979887\).
The [live project table](https://github.com/alexlegeartis/SpherePackings/blob/main/LOWER_BOUNDS_1_96.md)
lists the same occupied density \(2.26696032293109\cdot10^{-8}\).
The result does not concern other \(w\), partial fibers, arbitrary
non-dual translates, or all 47D packings. No density record changed.

The key new tool is to constrain **several projected norm shells at
once**, not to assume the first shell reaches every label. In
particular the tempting norm-18 three-sheet subgroup \(\{0,6,12\}\)
would have center density about 9719.18 if the label \(6\) lacked a
norm-6 lift. It does not: the small-neighbor argument below already
forces such a lift, independently of the large shell calculation.

## 1. Geometry and the small-neighbor argument

Primitivity and unimodularity imply \(x\mapsto x\cdot v\) on
\(\Lambda\) is onto \(\mathbb Z\); otherwise \(v\) is divisible in
\(\Lambda^*=\Lambda\). Slicing in unit steps gives
\(\operatorname{covol}(L)=\sqrt w\), since one unit of label has height
\(1/\sqrt w\). Projecting any \(x\in\Lambda\) gives an element of
\(L^*\), and its coset modulo \(L\) is determined by \(a=x\cdot v\)
modulo \(w\). Distinct residues are distinct classes: two projections
differing by an element of \(L\) have ambient difference, after
subtracting that element, in \(\Lambda\cap \mathbb Rv=\mathbb Zv\).
There are \(w\) classes, so these are all of \(L^*/L\). The horizontal
squared norm of a vector of ambient norm \(m\) and label \(a\) is
\(m-a^2/w\).

Suppose \(w=b^2t\), with \(t=2\) or 4 and \(b\) any positive integer.
Put \(K=\{x\in\Lambda:x\cdot v\equiv0\pmod b\}\) and
\(\Lambda'=K+\mathbb Z(v/b)\). The first index is \(b\) by
surjectivity, and the second is \(b\) by primitivity of \(v\) and
\(v\in K\). Pairings with \(v/b\) are integral on \(K\) and \(t\)
is even. Hence \(\Lambda'\) is even unimodular, also when \(b\)
is composite. Its section perpendicular to \(v/b\) is exactly the
union of the \(b\) canonical classes with labels multiples of \(w/b\).
Indeed the perpendicular projection of an original vector of label
\(j w/b\) subtracts exactly \(j v/b\), and the perpendicularity
condition on \(x+jv/b\) is precisely \(x\cdot v=-j w/b\).

Here is a general elementary harmonic-theta test. In an even
unimodular 48D lattice containing a vector \(r\) of squared norm
\(t=2\) or \(4\), there exists a nonzero vector perpendicular to
\(r\) of norm at most \(4\). Take
\(H(x)=(r\cdot x)^2-t\|x\|^2/48\).
The harmonic theta series is \(c\Delta E_4^2 E_6\) in weight 26,
with its \(q^2\) coefficient \(-48\) times its \(q\) coefficient.
If no root is perpendicular and roots exist, the \(q\) coefficient
is positive (each root contributes at least \(1-t/24>0\)).
If also no squared-norm-4 vector is perpendicular, the \(q^2\)
coefficient is nonnegative (each contributes at least \(1-t/12>0\)),
a contradiction. In the rootless case \(t=4\), the harmonic theta
series vanishes, whereas the shell containing \(r\) would have strictly
positive harmonic sum if no perpendicular vector existed.

The original \(L\) has no squared norm below 6, so this necessarily
supplies a **nonzero subgroup phase**. At \(w=18,b=3,t=2\), a vector
in either nonzero phase of norm at most 4 lifts to original norm at
most \(4+6^2/18=6\), by choosing the shortest normal label \(6\) or
\(-6\). This is exactly a norm-6 lift, so the projected three-sheet
minimum is exactly 4. At \(w=16,b=2,t=4\), the order-two phase has
squared clearance at most 4 by the same argument.

## 2. Synchronized harmonic equations for norms 6, 8 and 10

Write \(N_m(a)=\#\{x\in\Lambda:\|x\|^2=m,\ |x\cdot v|=a\}\).
Comparing \(x\) with \(\pm v\), which are of norm different from
these shell vectors, gives \(a\le(w+m-6)/2\); Cauchy–Schwarz
also gives \(a^2\le mw\). These are finite ranges.
When \(a> w/2\) at the boundary, \(x\mapsto v-x\) is a bijection
between positive-label shells \((m,a)\) and
\((m+w-2a,w-a)\). Negating doubles both sides when labels are nonzero.
This yields the boundary equations used below, with smaller shell
norms chosen from 6, 8 and 10.

The shell counts \(52\,416\,000,39\,007\,332\,000,
6\,609\,020\,221\,440\) are the coefficients at \(q^3,q^4,q^5\) of
\(E_4^6-1440\Delta E_4^3+125280\Delta^2\); the latter expression
is fixed by the constant term and vanishing through \(q^2\).
The three shells are 11-designs: their positive-even low degree
harmonic theta series vanish, since division by \(\Delta^3\) yields
negative weight. The exact absolute-label moment equations are
\[
\sum_a N_m(a)a^{2j}
=|\Lambda(m)| \prod_{i=0}^{j-1}\frac{mw(2i+1)}{48+2i},
\qquad j=0,\ldots,5.
\]

For even degree \(\ell\) use the zonal harmonic
\[
 H_\ell(x;v)=\sum_{j=0}^{\ell/2}
            c_j(mw)^j(x\cdot v)^{\ell-2j},\quad
c_0=1,\quad
c_j=-\frac{(\ell-2j+2)(\ell-2j+1)}
               {2j(48+2\ell-2j-2)}c_{j-1},
\]
where \(m=\|x\|^2\); the recurrence is the Euclidean Laplacian
identity. Its weight-\((24+\ell)\) theta series lies in
\(\Delta^3 M_{\ell-12}(\mathrm{SL}_2(\mathbb Z))\).
Thus for each \(\ell=12,14,\ldots,34\), the three shell sums of
\(H_\ell\) lie in the span of the three initial coefficients of
the basis
\[
 \Delta^{3+c}E_4^a E_6^b,\quad 4a+6b+12c=\ell-12,\quad b=0,1.
\]
The checker constructs every relevant annihilating row of this span.
For instance degree 14 has zero span, and in degree 12 the
ratio across the three shells is \(1,-72,2484\). All these
linear constraints are enforced simultaneously.

We relax counts to nonnegative rationals. The polytope is bounded
by the known total counts. Exact row reduction leaves at most four
parameters for the five normal lengths; evaluating its vertices
gives rigorous lower bounds on each shell bin. Positive lower bounds
give **upper** witnesses to the possible phase distance, even if
other shorter ambient lifts are also present. Examples: for \(w=18\)
the \(|a|=9\) bin has at least \(3\,688\,740\) norm-10 lifts; for \(w=20\)
the norm-10 \(|a|=9,10\) bins have lower bounds \(242\,003\,200/13\) and
\(299\,928\). Thus the apparent four-sheet win at \(w=20\) using
\(\{0,1,10,11\}\) does not survive.

## 3. Uniform phase-size calculation

The following are upper bounds on the squared clearance, multiplied
by \(w\), in order \(a=0,1,\ldots,w/2\). The zero entry refers to
the within-period minimum, which is 6 by the positive norm-6 zero
bin. Any value above this can be truncated to \(6w\).

| \(w\) | Clearance numerators | Maximum certified upper center density |
|---:|---|---:|
| 12 | 72, 71, 68, 63, 56, 71, 60 | 5712.700939 |
| 14 | 84, 83, 80, 75, 68, 59, 76, 63 | 5544.756951 |
| 16 | 96, 95, 92, 87, 80, 71, 92, 79, 64 | 5373.352601 |
| 18 | 108, 107, 104, 99, 92, 83, 72, 95, 80, 99 | 5207.135183 |
| 20 | 120, 119, 116, 111, 104, 95, 84, 111, 96, 119, 100 | 5049.564178 |

For each threshold \(B\) among the clearance numerators, take the
Cayley graph on \(\mathbb Z/w\) with edges permitting that threshold.
An actual motif of squared minimum \(\mu\) is a clique for the next
table threshold at least \(w\mu\). Fixing one vertex at zero by
translation removes no possibilities. If the clique cap is \(s\)
at \(B\), the squared center density is at most
\(\frac{s^2}{w}(B/(4w))^{47}\). The checker makes all comparisons with
\(\delta_*^2\) as integers. In all five rows the maximum upper bound
is at two adjacent phases; these do attain the indicated squared
minimum \(6-1/w\) (ambient minimum 6 and positive norm-6 label 1).

Run [verify.py](verify.py) with ordinary Python 3. The proof uses
exact rational polytope vertices and integer clique checks, not a
numerical LP solver or shell-sampling assumption. For the design and
weighted theta background see [Nebe's Venkov survey](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf)
and the [published 48D ambient constructions and catalogue references](https://arxiv.org/abs/1312.3780).

See also [the exact obstruction to two-neighbor edges between extremal 48D lattices](NO_EXTREMAL_2_NEIGHBORS.md).

[Fixed-period local rigidity of the two adjacent 12-norm sheets](LOCAL_TRANSLATION.md).
