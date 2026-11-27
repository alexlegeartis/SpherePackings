# Independently rotated D8 fibres: the sharp 9D bound still holds

Continue the normalization of [the main proof](README.md): each horizontal fibre is a full translate of an isometric copy of the unscaled root lattice \(D_8\), of squared minimum 2 and covolume 2. Allow **completely independent orthogonal orientations** in the same horizontal \(\mathbb R^8\), with arbitrary shifts and real heights, including repeated heights. Then at squared separation 2 the exact optimum over this whole class is still \(\Lambda_9\)'s center density \(1/(16\sqrt2)\). No common period, commensurability, or finite orientation palette is imposed. This addresses the rotation variant of the \(D_8\)-fibered packing problem mentioned by [Salmon, arXiv:2607.25254v1, after Theorem 13](https://arxiv.org/html/2607.25254v1#S3.SS2). It does not assert a sharp two-point LP certificate, or an unrestricted 9D packing bound.

## A gap between unlike orientations

Let \(A=D_8\) after choosing axes of the first fibre, and \(B=U D_8\) be the oriented lattice of a second one. If \(A\ne B\), the vertical separation \(h\) obeys, for any total squared separation \(\mu\),
\[
h^2\ \ge\ \mu-2+\frac{59}{216}.                         \tag{3}
\]
Here equality of orientations means equality of lattices as sets; changes by an automorphism do not count.

**At least two coordinates vary modulo integers.** Let \(G\) be the closure of \(B\bmod\mathbb Z^8\) in the coordinate torus. If all coordinates are constant it is trivial, so \(B\subseteq\mathbb Z^8\) is an index-two sublattice. Any index-two sublattice is the kernel of a nonzero binary linear form on integer coordinates. Since \(B\) has no unit vector, the form is nonzero on each unit vector and is exactly the parity check defining \(A\). Thus \(A=B\), excluded.

Nor can just one coordinate \(i\) vary. In that case the projection of \(B\) deleting coordinate \(i\) is a full-rank sublattice of \(\mathbb Z^7\), of some integer index \(k\). Its kernel in \(B\) is rank one on the \(i\)-th axis. Taking a kernel generator followed by lifts of a projection basis shows that the positive axial generator length is \(2/k\), because \(B\) has covolume 2. Inner products in the rotated lattice are integral, so pairing with that generator forces the \(i\)-th coordinate of each vector of \(B\) to be rational. Its square is an integer because the total squared norm and the other coordinate squares are integers; hence the coordinate itself is integral, contradiction.

**Each varying coordinate pays an average quadratic loss.** For the phase difference \(p\) between the fibres, all pairs satisfy the packing inequality, so for every \(v\in B\)
\[
q_A(p+v):=d(p+v,A)^2\ge\mu-h^2.
\]
With the main proof's notation set \(f(\delta)=\delta(1-\delta)\), where \(\delta_i\) measures distance to a half-integer and \(d=\min\delta_i\). For any selected set \(S\) of at most four coordinates, the remaining at least four contribute \(\sum f(\delta_i)\ge4f(d)\ge2d\). The parity cost is at most \(2d\), so
\[
\sum_{i\in S} f(\delta_i(p+v)) \le 2-q_A(p+v)\le 2-\mu+h^2.
\]
The left side is a continuous function on the coordinate torus and this inequality passes to \(p+G\).

The normalized Haar probability on \(G\) projects in each coordinate to the normalized measure on its image, which is either a group of order \(q\) of equally spaced points, or the whole circle. For a group of order \(q\), the average of \(\operatorname{dist}(x,\mathbb Z)^2\) over any translated \(q\)-grid is at most
\[
\frac{q^2+2}{12q^2}.
\]
Indeed it is convex between successive locations where one node crosses a half-integer, so its maximum is attained with a node at \(1/2\). For even \(q=2k\), the distances there are \(1/2,0\), and twice each \(j/q\), \(1\le j\le k-1\); for odd \(q=2k+1\), they are \(1/2\) and twice each \((j+1/2)/q\), \(0\le j\le k-1\). Summing gives the displayed formula in both cases. For the whole circle the average squared distance is \(1/12\). Since \(f(\delta)=1/4-\operatorname{dist}(x,\mathbb Z)^2\), a nontrivial coordinate pays mean loss at least \(1/8\), and it pays at least \(4/27\) if the order is at least three or infinite.

If at least three coordinates vary, choose three of them to obtain average cost at least \(3/8\). If exactly two vary, their orders cannot both be two: then every rotated lattice vector has at most two half-integral coordinates and all others integral, and the integrality of its squared norm forces the number of half-integral coordinates to be divisible by four. This would make both coordinates trivial. Thus exactly two varying coordinates cost at least \(1/8+4/27=59/216\). Both cases prove (3).

## Three arbitrary orientations

If two of the three fibres have the same orientation \(L\), every pair in the triple involves a fibre of orientation \(L\). The actual horizontal difference set contains the phase difference plus \(L\). The assumed packing separation therefore gives the same lower bound on \(d(\text{phase difference},L)^2\) as in the aligned proof. Its three-layer argument applies in one coordinate frame for all three pairs, even if one fibre is differently oriented.

If the three orientations are all distinct, apply (3) to the two successive gaps. At \(\mu=2\) their sum is at least \(2\sqrt{59/216}>1\). Thus every normal interval of length less than one again contains at most two fibres, without orientation constraints. At most two fibres per unit asymptotic height, horizontal covolume 2, and radius \(1/\sqrt2\) give center density at most \(1/(16\sqrt2)\), attained by the aligned fixture.

For the inherited-radius range \(t=\mu/2\in[93/100,1]\), any triple has span at least the smaller of \(\sqrt{3t-2}\) and \(2\sqrt{2t-2+59/216}\). The former is the same-orientation-pair argument; the latter is for three different orientations. Both squared quantities are at least \(t^9\) on the stated interval, by convexity of \(t^9\) minus the relevant affine function and the exact endpoint checks in the verifier. The bound at this interval of scales still does not exceed the \(\Lambda_9\) density.

For horizontal dimensions 9–11, the aligned argument also applies if there are at most two allowed orientations in the entire family (by the same pigeonhole principle in every triple). This rotation proof imposes no claim about three or more different orientations in dimensions 10–12.

