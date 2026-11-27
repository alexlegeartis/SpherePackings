# Literal 44D P48p embedding: a sharp dual-phase ceiling

For the particular determinant-605 normal in \(P_{48p}\) that [Sun and Wang realize in their source](https://github.com/wcgbg/sphere-packing/blob/a0bc7f51444e42c3cbae0d2a910039cd50013500/sphere_packing/core/chlww25_sections.cc) for the [Chen et al. 44D construction](https://arxiv.org/html/2505.02394v1#S3.SS7), **every union of full perpendicular sheets with dual phases has center density at most 509.6189826428645…**, the published value. This normal is not claimed universal across embeddings. The theorem does not restrain phases outside the dual group or metric changes.

The five normal vectors \(u_1,\ldots,u_4,-\sum_i u_i\) have squared norm 6; their mutual dot products \(-2\) form a pentagonal cycle, the other dot products being \(-1\). Equivalently the Gram of the supplied four vectors is

\[
G=\begin{pmatrix}6&-2&-1&-1\\-2&6&-2&-1\\-1&-2&6&-2\\-1&-1&-2&6\end{pmatrix},\qquad
55G^{-1}=\begin{pmatrix}14&8&7&6\\8&16&9&7\\7&9&16&8\\6&7&8&14\end{pmatrix}.
\]

The normal is primitive: a nontrivial saturation index would have to be 11 because \(\det G=605\), hence the saturated rank-four lattice would have Gram determinant 5 while retaining minimum at least 6. Even the elementary volume bound for disjoint 4-balls of radius \(\sqrt6/2\) makes this impossible. Thus the perpendicular period has covolume \(\sqrt{605}\) by unimodular gluing.

We explicitly supply one lattice vector per antipodal pair of dual phases, lifted at a **shortest** normal representative. In 547 of the 605 classes (including zero) the lift has squared norm 6; in the other 58 it has norm 8. The latter witnesses do **not** claim that a norm-6 lift fails to exist. They only give upper bounds on the true phase distances, which suffices. If \(q\) denotes 55 times the squared norm of that shortest normal representative, these square-distance upper bounds are \(6-q/55\) and \(8-q/55\), respectively.

The witness membership test reconstructs the Pless symmetry ternary code \([I_{24}\mid S]\) from its signed Paley matrix over \(\mathbb F_{23}\), in the same orientation as the [authors' lattice builder](https://github.com/wcgbg/sphere-packing/blob/a0bc7f51444e42c3cbae0d2a910039cd50013500/sphere_packing/lattices/p48pq.cc). Raw vectors are divided by \(\sqrt{12}\). The integer bank is \(2x\), \(x\bmod3\) in the code and \(\sum x_i\) even; the half-integer bank is \(3\mathbf1+2x\) with odd \(\sum x_i\). The opposite even neighbor has the norm-4 all-minus-one raw vector, so this choice is the external minimum-6 \(P_{48p}\). The packed witnesses are encoded as one letter per coordinate, with letter g giving zero and successive letters changing the coordinate by one.

Using those upper bounds jointly makes an optimistic Cayley graph for each radius. For integer \(d=55\mu\) among the candidate thresholds, an improvement would require a clique of more than \(5(314/d)^{22}\) phases. Fix one at zero and use exact greedy-color branch bounds to disprove such cliques in its neighborhood. For example \(d=316\) permits at most two phases in the optimistic graph, whereas \(d=314\) permits five, attaining Chen et al. The verifier uses the weaker density-derived clique limits where sufficient. All comparisons use integers, not a sampled clique search. Larger separations than the maximum mixed upper bound permit only a single translate, and the norm-6 zero-class witness bounds its minimum.

The five displayed dual basis phases including zero have all pairwise normal squared distances at most \(16/55\); hence they attain the lower bound with squared minimum \(314/55\). This recovers the external density, without a new record or authorship claim. Run the companion verify44.py with t44_witnesses.json to replay the finite certificate.

