# Γ72 covering and the 79D two-level normal-lattice gate

**No new density record is proved here.** This note gives (i) an exact reduction of full-fiber 73D improvement to a single covering-radius question, (ii) an all-radii sharpness theorem for integral-phase 73–75D fibers, and (iii) a concrete min-2 normal geometry that has 4.17% headroom in 79D, but provably requires nonlinear horizontal phase assignments. The neighboring 79D normal geometry is ruled out even for nonlinear assignments.

We use [Nebe's Γ72](https://arxiv.org/abs/1008.2862) in the [exact polarized block coordinates](../gamma72_depth6_20261122/), with minimum squared norm 8 and covolume 1. In three copies of the standard Leech lattice \(\Lambda\), the metric is half their direct sum. The Γ72 subgroup is characterized by equality of the three \(A\)-labels and a zero sum of \(B\)-labels in the polarization of \(\Lambda/2\Lambda\). Its integral phase group \(\Lambda^3/\Gamma_{72}\) has size \(2^{36}\); it is a proper subset of all real translation classes of Γ72.

## 73D: one real hole suffices, and is necessary at the inherited radius

Write \(d(x)\) for distance from a phase \(x\) to Γ72, and \(R\) for its covering radius. **At radius \(\sqrt2\), arbitrary complete translated Γ72 fibers in 73D can beat the current 73D density if and only if \(R^2>6\)**. The comparison is the [certified determinant-two stack](../gamma72_depth6_20261122/), center density \(2^{73/2}/\sqrt2=2^{36}\) and occupied fraction \(4.26462994945094\ldots\times10^{-14}\).

If \(R^2\le6\), all consecutive normal levels have spacing at least \(\sqrt2\) by their horizontal distances; hence this center density cannot be exceeded. Conversely, take a single phase \(b\) with \(d(b)^2>6\), and any established integral phase \(a\) of squared depth 6. Use a periodic four-layer horizontal motif
\[
0,\quad b,\quad a+b,\quad a
\]
and normal gaps \(g,h,g,h\), where \(h=\sqrt2\) and \(g<h\) is sufficiently close to \(h\). Adjacent horizontal differences alternate \(\pm b,\pm a\), so their squared distances plus the normal ones are at least 8. Across two gaps the differences are \(\pm(a+b),\pm(a-b)\); both have positive distance, since their being zero would force \(d(b)^2=6\). Taking \(g\) sufficiently close to \(h\) makes these pairs safe as well; all longer separations already clear in the normal coordinate. The density ratio to the incumbent is \(2h/(g+h)>1\). This uses a single new real hole, not a deep triangle of three holes.

There is a quantitative version. The [integral phase-graph certificate](../gamma72_phase_graph_20261126/) has four mutually deep phases including 0. Among the other three choose \(a\) outside both open balls around \(b,-b\) of radius \(\sqrt6/2\) in the horizontal torus: each such ball contains at most one of the three. Thus \(d(a\pm b)^2\ge3/2\). If \(D=d(b)^2>6\), it suffices to take
\[
g=\max\{\sqrt{8-D},\,\sqrt{13/2}-\sqrt2\}.
\]
Three or more gaps are safe even at the smallest allowed \(g\); e.g. \(2g+h\ge2\sqrt{13/2}-\sqrt2>\sqrt8\). A future literal depth-\(13/2\) witness would already yield approximately a 7.18% improvement by this formula. No such new hole is asserted.

### Where such a hole must be

At standard squared minimum 4 the covering radius of the Leech lattice is **\(\sqrt2\)**, not 2. The two polarized sublattices in Nebe's construction are each isometric to \(\sqrt2\,\Lambda\). For any real block triple \(x\), select the first approximation in the block minimizing \(d(x_i,\Lambda)\), then select the second in its prescribed \(A\)-class, then the third in its prescribed \((A,B)\)-class. Their unhalved squared errors are respectively at most \(d(x_i,\Lambda)^2,4,8\), by the covering radii of \(\Lambda,\sqrt2\Lambda,2\Lambda\). Hence
\[
 d(x,\Gamma_{72})^2\le 6+\tfrac12\min_i d(x_i,\Lambda)^2\le7.
\]
In particular **a deeper hole cannot have even one integral block**, and a hole of depth \(6+\epsilon\) requires each block to have squared distance at least \(2\epsilon\) from its own Leech lattice. This rules out moving only one or two blocks away from the integral phase stock. This is a bound, not an exact value of \(R\). See also the [polarized-lattice decoding work](https://perso.telecom-paristech.fr/ciblat/docs_research/journals/ieee_tit2022_decoding.pdf), which gives background on the covering-radius problem.

Every integral depth-six phase is in fact a strict local maximum of the torus distance, with an explicit neighborhood. In the phase graph's notation \((a,b,t)\), its two relative \(A\)-labels are independent. For each chosen block, set another block zero and prescribe the remaining two nonzero \(A\)-labels in the corresponding triangle. Every norm-four vector with the first prescribed nonzero \(A\)-label can be used: the other required \(B\)-label must give norm 8, since norms at most 6 would yield total depth less than 6. The 24 antipodal rays of the first norm-four shell are an orthogonal frame. Indeed both the sum and difference of distinct nonantipodal representatives lie in the polarized sublattice of minimum 8, so their inner product is zero. Signs in the two occupied blocks can be chosen independently without changing binary phases. For a displacement \(u\) of halved-metric length \(s\), a block has ordinary length at least \(\sqrt{2/3}\,s\); its 24-frame contains a contact direction whose inner product in the halved metric is at most \(-s/6\), with all other block contributions made nonpositive. Therefore
\[
 d(c+u,\Gamma_{72})^2\le6+s^2-s/3
\]
for each integral depth-six phase \(c\). No displacement of length strictly between 0 and \(1/3\) around any such phase can be the required deeper hole.

## Integral phases: all scalar metrics in 73–75D

Fix the sphere radius \(\sqrt2\) by uniform normalization and use complete translates of \(\sqrt t\,\Gamma_{72}\) for \(t\ge1\). Allow arbitrary positions and occupancies of integral phase fibers in the \(k\)-dimensional normal space, \(k=1,2,3\), not necessarily periodic or lattice-based. The center density divided by \(2^{(72+k)/2}\) is at most the normal point density divided by \(t^{36}\). The [existing phase-graph proof](../gamma72_phase_graph_20261126/) gives effective squared depth at most 6, maximum size 4 for a set of pairwise depth-six phases, and a necessary nonzero first projected label on each difference of squared depth at least 5. Hence a set of pairwise depth at least 5 has size at most \(2^{12}\). If \(P_k=1,2/\sqrt3,\sqrt2\) and \(V_k\) denote optimal normal point density at separation one and unit-ball volume, respectively, the bounds on the normalized center density are
\[
 \frac{P_k}{t^{36}(8-6t)^{k/2}},\quad
 \frac{4\,2^k}{V_k t^{36}(8-5t)^{k/2}},\quad
 \frac{4096\,2^k}{V_k t^{36}(8-4t)^{k/2}},\quad
 \frac{2^{36} P_k}{t^{36}8^{k/2}} .
\]
The first uses normal minimum, the next two integrate local counts in open balls of diameter \(\sqrt{8-jt}\), \(j=5,4\), and the last sums the normal packing bounds over phases separately. Apply these bounds successively over the four intervals with cut points \(1,5/4,3/2,79/40,\infty\). Each \(t^{72}(8-jt)^k\) has at most one interior maximum, so endpoint tests suffice on each bounded interval. Exact comparisons in the checker show these densities are below the inherited comparators \(P_k/2^{k/2}\), except equality at \(t=1\) from the first bound. The comparators are attained by the prior determinant \(2,3,4\) constructions. Coincident normal placements at large \(t\) are allowed with distinct labels and counted with multiplicity. This does not bound arbitrary real horizontal phases.

## 79D: normal geometry with headroom, but the phase assignment is not solved

Take a completely filled seven-dimensional normal lattice as follows. Let \(K=D_3\) in the three last normal coordinates, with basis \(k_0=(1,1,0),k_1=(1,0,1),k_2=(0,1,1)\). Choose an elementary index-four overlattice \(P\subset\tfrac12K\), and biject its three nonzero cosets modulo \(K\) with the three nonzero deep shifts of \(D_4^*/D_4\). Over each \(y\in P\) put one copy of the corresponding shifted \(D_4\), using squared metric multiplier \(r=3/2\) in those four coordinates. D4 deep shifts have squared norm depth 1 and the implied minimum normal squared distance is 2. The normal squared covolume is
\[
 D_0=(2r^2)^2(\operatorname{covol}(K)/4)^2=81/16.
\]
This is **only a normal geometry**, not a 79D witness until horizontal phases are supplied. The [79D incumbent MW80 antipode](../creative_sections_20261116/PROOF.md) has two translates, minimum \(3235/234\), period determinant \(41^{12}(234/41)\), so it corresponds at the Γ72 normalization to
\[
 D_*=\frac{2^{79}\,41^{11}\,234}{4(3235/936)^{79}}
 =5.49378276\ldots .
\]
If the required shifts existed on this normal geometry, the prospective gain would be \(\sqrt{D_*/D_0}-1=4.17\ldots\%\).

**Type 0**, \(P=\langle k_0/2,k_1/2,k_2\rangle\). Use D4 shifts \(e_1\) and \((1/2)^4\) for the first two tail generators. Take four ordinary D4 basis coefficients followed by the three stated tail coefficients, so each normal point is indexed by \(\mathbb Z^7\). This is the exact geometry used in the checker. The 60 vectors with norm below 3 all have norm 2, and represent 30 distinct nonzero coefficient parities; the 180 vectors below 4 represent 90 parities. A horizontal projected coloring \(p:\mathbb Z^7\to\mathbb F_2^{24}\) is forced to have projected differences in the preceding certificate's set \(R\) on the short edges, and nonzero projected differences on all edges of length squared below 4. If \(p\) is additive, its kernel is a binary subspace disjoint from the 90 parities. There are exactly 38 options: zero and 37 dimension-one subspaces (checking pairs rules out any plane).

For each quotient, identify the distinct images of the 30 short parities, and connect two when their sum is another such image. This graph must map injectively to the induced compatibility graph of \(R\), which has induced links that are matchings and maximum codegree four. The exact check finds an element of degree two *inside a link* in 37 cases and codegree exceeding four in the remaining case. Thus **no additive integral-phase coloring works**, even when D4 phase assignment is allowed to vary with position inside each normal layer. This does not rule out arbitrary nonlinear colorings.

A concrete sufficient nonlinear problem uses the 128 verified phases of the prior rank-law package. Find two proper four-colorings of the infinite normal short-edge graph (squared length 2). Use their color pair as the rank-law matrix \(U\in M_2(\mathbb F_2)\) by an identification in which \(\det U=q_-(p_1)+q_-(p_2)\), with \(q_-\) zero at 0 and one at each of the three nonzero vectors of \(\mathbb F_2^2\). Such an identification is obtained with the two perpendicular anisotropic planes spanned by \(I,\left[\begin{smallmatrix}0&1\\1&1\end{smallmatrix}\right]\) and by \(\left[\begin{smallmatrix}1&0\\1&1\end{smallmatrix}\right],\left[\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right]\). Then require the color pair to distinguish all edges of squared length below 4, and add a three-bit subcolor to distinguish all remaining pairs of squared distance below 8 that share the color pair. The saved literal 128 phases realize squared depth 6,5,4 when respectively the matrix difference is rank one, rank two, or zero with different subcolor. This would exactly meet the distance requirements. No such coloring is claimed.

**Type 1**, the comparison plane in \(K/2K\) spanned by bit masks 1 and 6, with bit 1 denoting \(k_0\). This choice has a five-point obstruction even for arbitrary nonlinear labels. In doubled normal block coordinates, the checker exhibits five points with nine pairwise distances squared below 3 and the remaining one below 4. Anchor one of the universal vertices to projected phase 0. Another universal vertex has three distinct neighbors in its link in \(R\), two of which are both adjacent to the third. This violates the matching-link property. All projections are distinct because all squared separations are below 4. This statement is for the literal defined metric and plane.

There is also a low-scale obstruction for the whole tail template: let \(1\le r<3/2\) and replace \(K\) by any three-lattice. Minimum normal squared norm 2 forces \(\min K\ge2\). To be competitive \(\det K=4D/r^4<22\), because \(D<D_*<11/2\). Since \(2P\) has index two in \(K\), its shortest vector has squared norm at most \((8\det K)^{1/3}<6\) by the three-dimensional Hermite bound. This is not in \(2K\) (its minimum is at least 8); it gives a tail vector of squared length \(<3/2\) with nonzero D4-hole shift. In the old D4 layer choose one vertex from each opposite pair of the Delaunay orthoplex around that deep hole. All four mutual normal distances squared are \(2r<3\); their distances squared to the displaced point are \(r+\|y\|^2<3\). This would need five mutually depth-six phases, impossible. Thus reducing the D4 squared scale below \(3/2\) cannot repair the candidate within that template.

**Replay:** run node verify.js. The short normal enumeration includes every vector of squared norm below 8: if \(h,y\) are integer doubled coordinates, the numerator is \(3\|h\|^2+2\|y\|^2\) over 8, so \(|h_i|\le4, |y_i|\le5\). It checks the comparator inequality, exact nominal normal geometry, both phase-graph obstructions and the scalar endpoints. The inherited \(\Gamma_{72}\) link and codegree certificate is linked above; no new high-dimensional sphere enumeration is assumed.

