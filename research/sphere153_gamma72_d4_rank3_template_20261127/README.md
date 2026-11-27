# A 79D obstruction for the full scaled-D4, index-four, rank-three normal template

**No new density record.** This closes a concrete construction with 4.17% prospective headroom. The obstruction applies to arbitrary nonlinear integral Γ72 phase choices. It replaces an additive-only exclusion at the nominal metric by a tiny obstruction in the actual Euclidean graph, and extends it to all rank-three tail metrics in this isotropic-D4 template.

## Normalization, template, and theorem

Use [Nebe's Γ72](https://arxiv.org/abs/1008.2862) horizontally, with squared minimum 8 and covolume 1, in its three polarized Leech block coordinates. An *integral phase* means a shift belonging to the integral Leech-block superlattice, modulo Γ72. The 79D packing must have separation at least \(\sqrt8\) (sphere radius \(\sqrt2\)). For the normal site lattice take a standard \(D_4=\{x\in\mathbb Z^4:\sum_i x_i\ \text{even}\}\) layer with squared scale \(r>0\), over each site \(y\) of a rank-three lattice \(P\); \(K\subseteq P\subseteq K/2\) is a rank-three lattice and \(P/K\) has order 4. Here \(K/2\) in the inclusion means the dilation \(\tfrac12K\), not a quotient. Assign its four cosets bijectively to the four classes of \(D_4^*/D_4\), fixing zero, and use the assigned four-dimensional lateral shift. The isotropic scale and the \(D_4\) orientation are common to all layers; the choice of Γ72 phase at each individual site is entirely free. An arbitrary bijection fixing zero on these binary planes is additive, so the *normal stock* is indeed a lattice. No additivity at all is assumed of the horizontal Γ72 phase map.

Write \(D\) for the **squared covolume of the normal lattice** and \(\det K\) for the squared covolume of the tail section \(K\). The normal point density is \((2r^2\,\operatorname{covol}(K)/4)^{-1}\), since \(D_4\) has unscaled covolume 2. Consequently
\[
      D=\frac{r^4\det K}{4},\qquad \delta_{79}=\frac{2^{79/2}}{\sqrt D}
\]
for full occupation. The [incumbent 79D MW80 two-translate antipode](../creative_sections_20261116/PROOF.md) has squared minimum \(3235/234\), two translates, and period determinant \(41^{11}\,234\). Its exact equivalent normal determinant is
\[
      D_*=\frac{2^{79}\,41^{11}\,234}{4(3235/936)^{79}}
          =5.49378276\ldots .
\]
For the nominal \(r=3/2,K=D_3\), the stock has \(D_0=81/16\), 4.17% potential gain if all pairs could be colored. The [previous additive calculation](../sphere151_gamma72_covering_and_79_tail_20261127/) had left the unrestricted nonlinear assignment open on type 0.

**Theorem.** Every fully occupied packing in the stated template has \(D\ge729/128=5.6953125>D_*\). Thus it cannot reach the incumbent. Moreover, for \(r\ge3/2\) and tail minimum squared at least 2, the conclusion holds even with arbitrary deleted normal sites or arbitrary deleted centers within one integral-phase horizontal fiber at each normal site: its density is at most \(2^{79/2}/\sqrt{729/128}\). The theorem does not exclude other lateral geometries, nonintegral phases, or different horizontal fibers.

## The inherited finite graph

Here are precisely the input facts from the [verified integral Γ72 phase graph](../gamma72_phase_graph_20261126/), especially its [codegree replay](../gamma72_phase_graph_20261126/verify_extra.py). Every integral phase depth squared is an integer at most 6. Projection to the pair of relative \(A\)-labels gives \(p\in\mathbb F_2^{24}\); a zero projection has squared depth at most 4. A depth-six difference projects into a set \(R\) of 5,850 nonzero labels. In its graph, \(a,b\in R\) are adjacent when \(a+b\in R\). Every induced link is a matching, and every two distinct vertices have at most **four** common neighbors.

The necessary \(R\)-conditions are geometric, not a fitted list. For distinct nonzero \(A\)-labels \(a,b,a+b\), use two occupied Leech blocks and put the third block zero. Let \(M(u)\) be the possible \(B\)-labels of norm-four classes at a nonzero \(A\)-label \(u\). If one occupied block has \(B\)-label \(x\in M(u)\), the other has prescribed label \(t+x\). If \(v\cdot(t+x)=1\), the latter Leech residue has norm six, giving a \(4+6\) representation and horizontal squared distance 5, not 6. Hence \(v\) must annihilate all differences of \(M(u)\), for every ordered pair \(u,v\) of distinct members of the triangle; the affine equations on the common \(t\) must also be consistent. The inherited shell check produces the 975 unordered triangles (5,850 ordered labels) meeting these necessary conditions. Its exact adjacency replay gives each neighbor link seven disjoint edges and maximum codegree four. The proof below only uses this necessary superset \(R\), so no sufficiency for actual horizontal depth is presumed.

Consequently, in a putative full-fiber packing, **normal squared distance below 3** forces a depth-six difference and membership of the projected difference in \(R\); **normal squared distance below 4** forces distinct projections. This converts the distance windows below into purely local graph contradictions.

## Two tiny exact local normal windows

Use nominal \(K=D_3\) to name the vectors, with \(k_0=(1,1,0), k_1=(1,0,1), k_2=(0,1,1)\). Coordinates in the display are **doubled** in both blocks, \((h,y)\); under a changed tail metric \(Q\), squared normal distance is \(r\|h-h'\|^2/4+Q(y-y')/4\). The four doubled \(D_4\) shifts in their standard order are \(0000,(2,0,0,0),(1,1,1,1),(3,1,1,1)\). Type 0 denotes the plane spanned by binary masks \(1,2\) of the \(k\)'s, and type 1 the plane spanned by \(1,6\). In listing shifts for type 1 order the classes \(0,1,6,7\).

**Type 0 (eight vertices):**
\[
\begin{array}{c|rrrr|rrr}
 &\multicolumn4c{h}&\multicolumn3c{y}\\
 O&0&0&0&0&0&0&0\\
 A&-2&0&0&0&-1&-1&0\\
 B&0&0&0&0&-2&0&-2\\
 C&0&0&0&0&-2&-2&0\\
 D_1&-1&-1&-1&-1&-1&0&-1\\
 D_2&-1&-1&1&1&-1&0&-1\\
 D_3&-1&1&-1&1&-1&0&-1\\
 D_4&-1&1&1&-1&-1&0&-1
\end{array}
\]
At the nominal metric, all seven points are at normal squared distance 2 from \(O\). Also \(A,B\) each have squared distance 2 from every one of \(C,D_1,\ldots,D_4\). The pair \(A,B\) has squared distance 3, and the five points \(C,D_1,\ldots,D_4\) are mutually at squared distance 3. In particular the latter five project to **five distinct common neighbors** in \(R\) of two distinct vertices \(p(A)-p(O),p(B)-p(O)\). This contradicts codegree four. There is no need to classify 4-colorings or use parity-quotient fictitious edges. The exact same contradiction uses strict bounds below 3 for the edges just called short, and strict bounds below 4 for the edges of squared length 3. In the nominal metric it is stable even when the actual horizontal shifts differ from integral shifts (after a common translation) by lengths \(<(\sqrt6-\sqrt5)/2\): nondeep integral pairs cannot be raised from at most \(\sqrt5\) to the required \(\sqrt6\), nor zero-projection pairs from at most 2 to the required \(\sqrt5\).

**Type 1 (five vertices):** the order-\(1,6\) witness from the [previous note](../sphere151_gamma72_covering_and_79_tail_20261127/) is
\[
 (0,0),\quad(0,(-2,-2,0)),\quad(0,(-2,0,-2)),
\quad((-1,-1,-1,1),(-2,0,0)),\quad((-1,-1,1,-1),(-2,0,0)).
\]
The first zero within each pair is in the four-dimensional block. All nine pairs except the last pair must have squared distance below 3, and that last pair below 4. Anchor the projection at the first point. The other four project to four distinct vertices of \(R\), with every edge present except possibly the last pair. The link at the second projected point contains three vertices of which one is adjacent to both others. This is not a matching.

The signed coordinate permutations of the **nominal** \(D_3\) have exactly two orbits on the seven planes, of sizes four and three, represented by type 0 and type 1. The checker computes them by literal action on \(D_3/2D_3\). The permutations are only a convenient way to transform the *vectors* of each witness to an arbitrary plane; they are not assumed to preserve an actual tail metric. The three nonzero \(D_4^*/D_4\) classes can be permuted independently by isometries of \(D_4\): change one coordinate's sign to exchange the two half classes, and use the orthogonal four-by-four Sylvester Hadamard divided by 2 to send the odd integer class to the all-halves class. Thus any actual assignment of lateral class names is allowed.

## Every competitive rank-three tail preserves the windows

First recall the low-scale branch, which is also in the previous geometry note but holds for any rank-three tail. Full occupation forces \(r\ge1\) since same-layer \(D_4\) roots have squared normal norm \(2r\), and forces \(\min K\ge2\) in squared norm since distinct normal sites can have horizontal depth at most 6. Suppose \(D<729/128\). For \(1\le r<3/2\), \(\det K=4D/r^4<729/32\). The lattice \(2P\) has index two in \(K\), so by the exact three-dimensional Hermite constant \(2^{1/3}\) it has a nonzero vector of squared norm at most \((8\det K)^{1/3}<6\). This vector is outside \(2K\); half of it is a nonzero tail shift class represented by a vector of squared length \(<3/2\). In the original \(D_4\) layer, choose one vertex from each opposite pair of the Delaunay orthoplex around the associated deep hole. These four vertices are mutually at squared normal distance \(2r<3\) and at squared distances \(<r+3/2<3\) from the shifted point. Five mutually deep projections are impossible by the matching-link property.

For \(r\ge3/2\), the same Hermite constant and \(\min K\ge2\) give \(\det K\ge4\), so
\[
r^4\le D<729/128,\qquad
r<31/20,\qquad
\det K=\frac{4D}{r^4}<\frac92.
\]
We need only a weak, quantitative reduction of the metric to the FCC shape. It is stated next, with a finite exact certificate; it does not assume \(K\) itself is isometric to \(D_3\).

### Rank-three conorm lemma

If \(K\) has minimum squared at least 2 and squared covolume \(<9/2\), there is a basis identification of its underlying \(\mathbb Z\)-module with \(D_3\) such that every nominal vector of squared Euclidean norm \(2,4,6\) has actual squared norm, respectively, **less than \(3,29/5,49/5\)**.

Start with any superbase \(v_0+v_1+v_2+v_3=0\), where any three form a basis. If \(v_i\cdot v_j>0\), replace \(v_i\) by \(-v_i\), keep \(v_j\), and replace each of the other two by \(v_k+v_i\). This is again a superbase, via unimodular additions in the three-vector basis omitting \(j\), while the sum of squared lengths drops by \(2v_i\cdot v_j>0\). There are only finitely many bounded lattice vectors, so the procedure reaches an obtuse superbase. Put \(p_{ij}=-v_i\cdot v_j\ge0\), and choose the pair of opposite edges with smallest \(p_{ij}+p_{kl}\) to name the diagonals \(02,13\). This is a relabeling, not an optimization assumption. The conorm order is
\[
p=(p_{01},p_{12},p_{23},p_{30},p_{02},p_{13}).
\]
We have \(p_{02}+p_{13}\le p_{01}+p_{23},p_{12}+p_{30}\) (both inequalities). Also every nontrivial cut of the four vertices has weight at least 2, since that weight is the squared norm of the corresponding sum of superbase vectors. For integer coefficients \(a_i\), the exact squared norm is \(\sum p_{ij}(a_i-a_j)^2\). Identify the \(v_i\) in order with the nominal FCC superbase \((1,1,0),(-1,0,1),(1,-1,0),(-1,0,-1)\). Every three of the latter form a basis of the standard \(D_3\), so this is an integral isomorphism. A nominal \((x,y,z)\in D_3\) has superbase coefficients \(((x+y+z)/2,z,(x-y+z)/2,0)\).

For each nominal vector of the three lengths, add to these linear constraints the closed condition that its actual squared norm is at least \(3,29/5,49/5\), respectively. The resulting conorm polyhedron has nonnegative coordinates, and the reduced Gram of \(v_0,v_1,v_2\) is positive definite by the cut constraints (the positive-weight graph is connected). On each such polyhedron the minimum determinant is found among its vertices. For a direct justification, at a point whose active constraints have rank less than six, move along a null direction of these constraints. At least one endpoint is finite since all coordinates are nonnegative. If both are finite, one endpoint has no larger determinant by concavity of log determinant on positive definite Grams. If only one is finite, the direction from that endpoint toward infinity is coordinatewise nonnegative, so the Gram increases in positive-semidefinite order and the finite endpoint again has no larger determinant. The new active constraint increases rank. Iterating proves the claim without ignoring recession directions.

This is a small rational calculation: for each of the 6, 3, and 12 nominal vectors up to signs, enumerate six independent equalities from the 16 constraints (6 signs, 7 cuts, 2 opposite-pair inequalities, 1 norm cutoff); the checker uses exact fraction-free elimination and verifies feasibility. The minimum determinants over each respective family are
\[
\frac92,\qquad \frac{779}{150},\qquad\frac{1053}{200}.
\]
They are all at least \(9/2\), contradicting the assumed strict volume bound. For example the first value is attained for one nominal root with conorms \((1/2,1,1/2,1,1/2,1/2)\). This proves the conorm lemma.

### Substitution into the two windows

The literal vectors above, after their nominal signed permutations, still belong to the actual normal stock by their integer labels. All short type-0 edges use either double a nominal tail root and no \(D_4\) change, or a nominal tail root and a doubled \(D_4\) displacement of squared norm 4. Their squared normal distances are therefore \(<3\) and \(<31/20+3/4<3\), respectively. The other edges needed for distinctness use either no tail change and \(D_4\) displacement of doubled squared norm 8, or a nominal tail norm-six vector and doubled \(D_4\) displacement of squared norm 4. Their distances are \(<2(31/20)<4\) and \(<31/20+(49/5)/4=4\), respectively. In type 1 the short edges either use double a nominal tail root with no \(D_4\) change or use a nominal norm-four tail vector and doubled \(D_4\) norm 4; the latter are \(<31/20+(29/5)/4=3\). The exceptional pair is \(<2(31/20)<4\). The two local graph contradictions apply throughout the volume range.

For arbitrary deletion in the \(r\ge3/2,\min K\ge2\) branch, every translate of the eight- or five-point forbidden window must omit at least one full site. This stays true in *weighted occupancy* if only some horizontal points of the fiber at a site are kept: if no conflict existed between any two whole phases of the window, it would contradict the same local graph. Therefore at least one pair of phases has a forbidden difference, which gives a bijection between the two horizontal lattices pairing centers too close to retain both. Their combined asymptotic count in a common large horizontal cube is at most that of one fiber, up to a boundary error uniform in the shifts. Thus each window has total weight at most 7 or 4, respectively. Averaging translates bounds relative occupation by \(7/8\) or \(4/5\), valid also without periodicity. At \(D<729/128\), because \(D\ge r^4\ge81/16\), the effective determinant after deletion is at least
\[
       D/(7/8)^2\ge324/49>729/128.
\]
If the original \(D\) is already \(\ge729/128\), deleting cannot increase density. This completes the stated lower bound for the normal determinant equivalent.

Run `node verify.js` for the two literal windows, plane orbits, exact conorm polyhedra and comparator arithmetic. The published Γ72 graph replay is a separate inherited certificate. No native enumeration of Γ72 lattice vectors beyond that prior certificate is needed.
