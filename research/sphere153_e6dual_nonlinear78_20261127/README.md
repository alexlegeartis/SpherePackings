# The 78D scaled dual-\(E_6\) carrier: nonlinear full-stock and period-two obstructions

**No new record.** This closes the full-stock nonlinear loophole left in the [earlier scaled dual-\(E_6\) proposal](../gamma72_phase_graph_20261126/), and exactly excludes even nonlinear deleted support of period two in its projected colors. It is separate from the [79D rank-three tail theorem](../sphere153_gamma72_d4_rank3_template_20261127/).

Use complete translates of [Nebe's Γ72](https://arxiv.org/abs/1008.2862), covolume 1 and squared minimum 8, at sites of the six-dimensional normal lattice \(L\) with Gram \(G=(3/2)C^{-1}=H/2\). Here \(C\) is the \(E_6\) Cartan matrix on vertices \(0,\ldots,5\), with \(-1\) edges \(01,12,23,34,25\). This is the dual \(E_6\) normal lattice scaled to squared minimum 2. The horizontal shifts are integral polarized-block phases, with no additivity assumed. The required full packing separation squared is 8. The exact competitor is the [current 78D half-root Γ72–D4 lift](../gamma72_root_tilt_20261122/), whose normal determinant is \(27/4\); see [the project lower-bound table](../../LOWER_BOUNDS_1_96.md). Our stock has determinant \(243/64\). Its full-occupancy nominal center density \(2^{39}/\sqrt{243/64}\) is \(4/3\) times the competitor, so deleting down to relative occupation \(3/4\) would exhaust its headroom.

The two exact statements are:
* On the full nominal \(L\), no arbitrary nonlinear integral phase assignment works. For an arbitrary deleted support the **necessary** relative occupation bound from this proof is \(7/8\); it is **not** small enough to rule out a gain.
* If **both** the support and the projection of the integral phase map to the two relative \(A\)-labels are periodic under \(2L\), their relative occupation is at most \(44/64=11/16\). There is no requirement that the final \(B\)-label also have period \(2L\). Center density is at most \(11/12\) of the comparator. This does not settle longer-period or aperiodic thinnings or continuous phases.

## A real eight-point witness

We use the [published Γ72 phase graph](../gamma72_phase_graph_20261126/) and its [codegree replay](../gamma72_phase_graph_20261126/verify_extra.py). If two normal sites are closer than squared distance 3, their projected phase difference must belong to the set \(R_\Gamma\) of projected deep labels. If their squared distance is below 4, their projected labels must be distinct. In the graph on \(R_\Gamma\) whose edges are differences in \(R_\Gamma\), two distinct vertices have no more than four common neighbors. These are necessary conditions for true phases; no graph-to-phase sufficiency is asserted.

Take the following integer coefficient vectors for \(L\):
\[
\begin{array}{c|rrrrrr}
 O&0&0&0&0&0&0\\
 A&1&0&0&0&0&0\\
 B&1&0&0&0&0&-1\\
 D_1&0&1&0&0&0&-1\\
 D_2&1&-1&1&0&0&-1\\
 D_3&1&0&-1&1&0&0\\
 D_4&1&0&0&-1&1&0\\
 D_5&1&0&0&0&-1&0 .
\end{array}
\]
All seven nonzero vertices have squared normal distance 2 from \(O\); \(A,B\) both have squared distance 2 from each \(D_i\). The pair \(A,B\) and all pairs among the \(D_i\)'s have squared distance 3. Anchor at \(O\). The five \(D_i\)'s have distinct projected labels, and these are five common neighbors of two distinct vertices of \(R_\Gamma\), a contradiction. Thus every translate of the physical window occupies at most seven sites; averaging gives the relative occupation cap \(7/8\) even without periodicity.

## Period two turns every parity into a usable constraint

In \(L/2L=\mathbb F_2^6\), all 63 nonzero parities have shortest squared length either 2 or 3. Indeed \(H^{-1}=C/3\) has diagonal entries \(2/3\), so any vector with \(H[z]\le6\) has \(|z_i|\le2\), and the checker enumerates that exact box. There are 27 short parities, those with \(q(z)=H[z]/2\bmod 2=0\) apart from zero, and 36 medium parities with \(q=1\). The form is nondegenerate since \(H C=3I\). Write \(\mathcal R\) for the 27 short normal parities. When both normal occupancy and the projected \(A\)-labels have period \(2L\), every short parity difference forces \(R_\Gamma\)-membership (use a true shortest pair), and all pairs of distinct occupied parities have distinct projected labels by a short or medium shortest pair.

Make an undirected graph on all 64 parities using displacement set \(\mathcal R\). At each occupied anchor, its occupied neighbors must avoid every 7-set formed by two neighbors and their five common neighbors *within the anchor's neighborhood*. Each such 7-set together with the anchor is a forbidden 8-window, by the same codegree contradiction. The [checker](verify.js) finds 216 such neighbor sets, and 13,824 distinct 8-windows under all anchor translations.

Why can one focus on an anchor missing at most nine of its neighbors? The full parity graph has degree 27 and least eigenvalue \(-5\): for a nontrivial binary additive character \(\chi\), the eigenvalue is
\[
 \sum_{z\in\mathcal R}\chi(z)=\tfrac12\sum_z (-1)^{q(z)}\chi(z)-1=\pm4-1.
\]
The last sum has absolute value 8 by nondegeneracy: squaring and translating one summation variable leaves only the zero difference by character orthogonality. For an occupied set of size \(M\) with \(e\) induced edges, decompose its indicator into the constant component and its orthogonal complement. Then \(2e\ge M^2/2-5M\), so the average number of absent root neighbors of an occupied vertex is at most \(32-M/2\). If \(M\ge45\) (at most 19 holes), at least one occupied vertex has at most nine holes among its 27 neighbors. Translate it to zero.

## Exact finite local classification, then forced outside holes

Enumerate the deletion subsets of the 27 neighbors of sizes 5 through 9 that meet all 216 neighbor windows. There are none of size 5, hence none smaller (a smaller one could be enlarged to size 5). Use the 36 orthogonal transvections \(x\mapsto x+B(x,a)a\), where \(q(a)=1\) and \(B\) is the polar bilinear form. These fix the anchor and preserve the parity graph, hence all the window conditions. There is no assumption that they generate every automorphism; the calculation explicitly covers all subsets by orbits under just these permutations.

For a representative with local deletion set \(D\), all other neighbors and the anchor are occupied. Any global window avoiding \(D\) and containing just one vertex outside the closed anchor-neighborhood forces deletion of that single outside vertex. Collect all such forced vertices. Among windows avoiding both the local deletions and these forced outside vertices, take a disjoint subfamily of their restrictions to the outside. These pairwise disjoint required hit sets cost one additional deletion each, independently of all already counted holes. The certificate keeps an explicit greedy disjoint subfamily (no optimization assumption); a zero restriction would already be impossible. The exact aggregate replay is:

| Missing neighbors | Local hitting sets | Transvection orbits | Minimum certified global holes |
| ---: | ---: | ---: | ---: |
| 5 | 0 | 0 | n/a |
| 6 | 72 | 1 | 26 |
| 7 | 1,728 | 3 | 26 |
| 8 | 19,440 | 7 | 24 |
| 9 | 136,350 | 22 | 22 |

Each viable local case forces more than 19 holes, contradicting the assumed \(M\ge45\). Therefore \(M\le44\), giving the stated density ratio \((44/64)(4/3)=11/12\).

Run `node verify.js`. The code checks the exact shell box, literal physical witness, the finite neighborhood and orbit coverage, the forced/independent deletion obligations and density arithmetic. The analytic step is the one-character eigenvalue calculation above. The inherited Γ72 graph replay is separate.
