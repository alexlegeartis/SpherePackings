# Sakamoto 85D projections: exact obstruction in 52 low-determinant charts

In **each** of Sakamoto's two rank-88 extremal even unimodular lattices
\(L,L_0\), choose one explicit normal triple for each of the **26**
possible Gram types of independent minimum vectors with determinant
at most 384. All 52 resulting 85-dimensional charts have center
density **at most \(4/5\) of the BN85 benchmark** for every finite
motif of projected phases. This extends the earlier two selected
\(2A_3\) charts to 52 selected charts. It does not classify *embeddings*
of each Gram type, forbid other projections, or improve any packing
density bound.

The parent Gram matrices and the 44 seeds in each lattice come from
[Sakamoto's pinned exact reproducibility data](https://github.com/RyotaroSakamoto0113/rank-88-extremal-lattice-computations/tree/628ce2a7ac23acf8d4c295677e4310c9ec5b21d9).
The parent minimum 8 and unimodularity are by
[Sakamoto's theorem](https://arxiv.org/abs/2609.15633v2).
The standalone phase-graph replay checks every Gram type and all
trial norms as exact integers; it does not assume shell completeness.

## Scope of normal charts

Let \(K=\sum_i\mathbb Z v_i\) for independent norm-8 normals
\(v_1,v_2,v_3\), with Gram \(G_K\) and determinant \(D\).
Every off-diagonal entry is an integer of magnitude at most 4
(apply the minimum to \(v_i\pm v_j\)). By sign and permutation
sort their magnitudes \(a\le b\le c\) and put any negative sign on
the last entry. If all three are nonzero the sign of their product
is the invariant; with a zero only the nonnegative case is needed.
Test the norm bound on \(v_1\pm v_2\pm v_3\) and
\(0<D=512+2abc-8(a^2+b^2+c^2)\le384\) (taking \(c\) signed).
The replay obtains exactly 26 necessary Gram types, all realized in
both actual parents by the supplied explicit triples. An embedding
of any such \(K\) is primitive: an index-at-least-two saturation
would have determinant at most 96 and minimum at least 8, whereas
the exact rank-three Hermite constant gives at most
\(2^{1/3}96^{1/3}<8\).

The horizontal \(T=L\cap K^\perp\) has covolume \(\sqrt D\)
by primitivity and unimodularity. The group of projected phases
modulo \(T\) is \(K^*/K\) of size \(D\).
For a parent vector \(x\), its dot-label vector \(b=(v_i\cdot x)_i\)
specifies its phase by \(\operatorname{adj}(G_K)b\) modulo \(D\).
All group elements are obtained by the three adjugate generators.
Its projected squared norm is
\[
 \|P_{K^\perp}x\|^2=\|x\|^2-b^\mathsf{T}G_K^{-1}b=8-h/D,\quad
 h=D(8-\|x\|^2)+b^\mathsf{T}\operatorname{adj}(G_K)b.
\]
Thus a literal vector gives a blocker at that score for both the
phase and its negative. No assumption about its being shortest
in the phase is needed. Each chart includes an actual parent
norm-8 row with *zero dot-label*, placing a norm-8 vector in
its horizontal period; \(m>8\) is therefore impossible.

## Exact motif upper bound

For \(m=8-h/D\ge0\), two phase sites with difference having a
checked blocker of score larger than \(h\) cannot coexist.
The compatibility graph is invariant under phase translations,
so maximum clique can be searched through the zero phase.
The replay uses exact bitsets and greedy proper coloring
as branch bounds. Let \(M\) be the graph's maximum clique.
The ratio of squared center density to squared BN85 density
is at most
\[
 \frac{4M^2(8D-h)^{85}}{D(8D)^{85}}.
\]
Here the benchmark is \(2^{41.5}\), the threshold previously used
for these charts. Between successive blocker scores the graph
does not change and this ratio decreases. When \(M=D\) already
bounds it below \(16/25\), all larger \(h\) are safe without a graph
search. Every integer score checked in all 52 charts passes
the exact \(16/25\) comparison, which proves the \(4/5\) cap.
The score list is not a claim about possible actual minima.

## Regenerate short vectors and replay

Run OPENBLAS_NUM_THREADS=1 python -B verify.py (NumPy, one CPU).
The input compresses the two Grams, the 44 seeds and two isometries
per parent, and the 52 literal normal triples. It does **not** ship
the large search banks.

The replay starts from the exact order-\(5\)/order-\(23\) orbit
bank of the seeds and grows it by \(u-v\) between two norm-8
vectors with dot 4. The resulting norm-8 rows are expanded under
the checked isometries. Differences at dot 3 give norm-10 blockers.
Both norm lists are independently checked by exact integer
quadratic forms; the sampled pairs and seeds are only a method to
find enough literal vectors. Matrix products include explicit
integer overflow guards. The floating dot proposals use operands
and all partial sums small enough to be exactly represented, and
the resulting rows are independently certified by their integer
norms. No floating value is used in the density comparison,
in the phase graph, or to certify absence of short vectors.

The proof gives no upper bound for other embeddings of these
26 types, higher determinant triples, normals of larger norm, or
general 85-dimensional packings.
