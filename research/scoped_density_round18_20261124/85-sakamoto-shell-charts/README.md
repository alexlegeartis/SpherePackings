# Shell-aware obstruction in two Sakamoto 85D projection charts

For each of the two [Sakamoto even unimodular rank-88 lattices](https://arxiv.org/abs/2609.15633v2)
one explicit 3-dimensional normal space is closed: **no choice of
projected phases** in that space improves the 85D center density
\(2^{41.5}\). This is a bound on actual projected distances, not
the uniform-loss certificate. Other charts remain open.

The exact parent Grams come from data/L/lattice.json and
data/L0/lattice.json at [the pinned source revision](https://github.com/RyotaroSakamoto0113/rank-88-extremal-lattice-computations/tree/628ce2a7ac23acf8d4c295677e4310c9ec5b21d9/data).
Parent minimum and unimodularity are inherited from Sakamoto.
Each three-vector normal Gram is \(G_3=4(I+J)\) of determinant 256.
Its integral span \(K\) is primitive: a proper saturation would
have Gram determinant at most 64 and squared minimum at least 8,
contrary to the rank-three Hermite bound. The orthogonal period
has covolume 16 and squared minimum 8, the latter by the
[26-chart harmonic-shell certificate](../83-87-antipode/SHELL.md).

The phase group is \(K^*/K\), of order 256: normal labels differing
by an element of \(K\) give the same horizontal phase, and
conversely. Represent a dot-label \(b\in\mathbb Z^3\) by
\((4b_i-\sum_j b_j)_i\) modulo 16, since
\(G_3^{-1}=(4I-J)/16\).
A literal parent vector of squared norm \(e\in\{8,10\}\) and
normal label \(b\) gives horizontal squared distance
\[
 e-\bigl(4\sum_i b_i^2-(\sum_i b_i)^2\bigr)/16=8-h/16.
\]
Its residue and negative residue are forbidden as phase
differences at any minimum larger than \(8-h/16\).
The inputs have positive-loss blockers on 132 and 192 residues,
respectively. No shell exhaustion is assumed.

At proposed squared minimum \(8-h/16\), the motif is a clique
in the Cayley graph keeping residues whose checked losses
are at most \(h\). A largest clique may be translated to include
zero. Exact bitset branching with proper greedy coloring gives
these caps at the relevant thresholds:

| \(h\) | \(L\) | \(L_0\) |
|---:|---:|---:|
| 0 | 8 | 5 |
| 3 | 11 | 6 |
| 4 | 13 | 8 |
| 8 | 16 | 12 |

The ratio of squared center density to squared target is at most
\(M^2(128-h)^{85}/(64\cdot128^{85})\), where \(M\) is the cap.
At fixed graph the bound decreases with \(h\); once even 256
phases cannot beat the target, further thresholds need no search.

For a future chart, let \(q_c\) be the least squared normal
norm in phase class \(c\). Replacing a normal label by its
nearest representative does not change the horizontal projection.
Thus its horizontal depth is exactly \(8-q_c\) when the norm-8
shell meets \(c\), and is at least \(10-q_c\) otherwise.
This gives a constructive clique search when shell occupancy
is complete. The present one-sided witnesses cannot certify
other residues empty.

Replay: python -B verify_actual.py (standard library only).
