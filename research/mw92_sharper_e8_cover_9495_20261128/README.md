# A sharper certified rank-eight cover improves the 94D and 95D MW92 lattices

| Dimension | Occupied density at least | Strict gain against previous repository proof |
|---:|---:|---:|
| 94 | 8.08605570564669e−21 | 0.00927% |
| 95 | 2.13790360855475e−21 | 0.01058% |

Both results are lattice lower bounds, with attained squared minimum 16. Their horizontal shears are existential, not literal 94D and 95D basis matrices. These small improvements are over the exact previous [94D](../mw92_e8_tube94_20261128/README.md) and [95D](../mw92_e8_tube94_20261128/95D.md) proofs, not claims of optimality.

## Certified rational eight-dimensional cover

The inspiration is [Dutour Sikirić, Schürmann and Vallentin, Duke Math. J. 142 (2008), Table 2](https://arxiv.org/pdf/math/0601084v4), where their approximate covering density is 3.142202. The decimal Gram file currentbest8.txt in the [Schürmann–Vallentin arXiv source archive](https://arxiv.org/src/math/0405441v4) (SHA-256 f7e4f278526ef709c2d802a64fe1d26df20e70dd3ff98a8273cb85b2e8a3d5bd) gives a near permutation-symmetric physical metric when converted from Magma's standard E8 basis. We use the following **rational** metric on the ordinary E8 lattice and certify it independently; the literature decimal is not an assumption.

Set \(H=K/10^9\), with \(K_{00}=797943325\), \(K_{0i}=K_{i0}=1030483\) for \(i>0\), \(K_{ii}=1120006184\) for \(i>0\), and \(K_{ij}=-52646518\) for distinct \(i,j>0\). Work in physical E8 coordinates, so E8 is the union of integer and half-integer vectors of even coordinate sum and has ordinary covolume 1. For \(A=K_{00},B=K_{11},C=K_{01},D=K_{12}\) set \(U=B-D\), \(V=B+6D\), \(T=AV-7C^2\). On the six-dimensional tail sum-zero space K acts by U, and in the orthogonal span of the first coordinate and normalized tail sum its matrix is \(\begin{pmatrix}A&\sqrt7 C\\\sqrt7 C&V\end{pmatrix}\). Thus \(A,U,T>0\) verify positive definiteness and \(\det H=U^6T/10^{72}\). The transformed Euclidean lattice \(H^{1/2}E_8\) has covolume \(\sqrt{\det H}\).

Here is the entire radius check. For a point x in the H-metric Voronoi cell of 0, put \(y=4\cdot10^9 Hx\); every doubled lattice vector v gives \(y\cdot v\le v^\mathsf T K v\). Tail-coordinate permutations preserve both K and E8, so assume \(y_1\le\cdots\le y_7\). Inequalities from \(v=\pm4e_i\) give \(|y_0|\le4A\) and \(-4B\le y_1\le\cdots\le y_7\le4B\), a product of an interval and a seven-dimensional simplex with the 16 initial vertices in the checker. Intersect with the 244 inequalities of the nonzero vectors v whose coordinates have common parity, sum divisible by four, ordinary squared norm at most 32, and sorted tail. Using just these valid inequalities already gives an outer enclosure of the true chamber. Rearrangement also shows these are the strongest instances of their respective tail-permutation orbits.

Each rational clip keeps old vertices on the retained side and intersects the cutting plane with segments between strictly separated old vertices having at least seven common tight constraints. This cannot miss a true edge in dimension eight: its containing facets include seven linearly independent active normals, all drawn from the recorded inequalities. Any additional tested segment has endpoints in the old convex region and its intersection cannot enlarge the clipped region. Inductively all true vertices survive, although some retained points need not be vertices. Points created on an intersecting segment inherit the common active constraints plus the cutting plane (both endpoints satisfy every old inequality); repeated coordinates are merged by unioning active sets. There are 806 retained points. The outer polytope is bounded by the initial chamber; it is full-dimensional, since sufficiently small strict ordered coordinates satisfy every non-chamber bound strictly. Convexity implies \(x^\mathsf T Hx=y^\mathsf TK^{-1}y/(16\cdot10^9)\) is bounded by its maximum on these retained points. The exact maximum is
\[
 R^2=\frac{4012913824940349978129273030643634963}{4012913856863815325464678688000000000}<1.
\]
With \(\pi<3.14159265359\), the verifier proves \(\pi^4(R^2)^4/(24\sqrt{\det H})<3.143\) by squaring and integer comparison. No Voronoi software or published floating output is relied upon.

## Transfer to the two packing bounds

The [previous proof](../mw92_e8_tube94_20261128/README.md) constructs a literal sublattice \(S\cong\sqrt8E_8\) in MW92 and an ordinary periodic equal-ball replacement of its weak tubes; only that independently verified embedding is inherited. Set \(p=14/3\). Rescale the certified cover to squared radius \(2p\); its weak-ball replacement count per original point is less than
\[
 \frac{3.143\,8^4}{(\pi^4/24)(2p)^4}<209/500
\]
using \(\pi>3.14159\). As before, because the count inequality is strict, first shrink slightly to give radius slack, then approximate a basis by rational combinations in the basis of S. The operator norm of the basis change tends to one, so both inequalities persist. The resulting covering is commensurable with S; replacing each weak S-tube and transferring to a common cell inside MW92 gives the stated point-count rate. This preserves the exact counting and whole-small-ball coverage of the predecessor.

Write \(A_0=11.00098526\) for the inherited MW92 Coxeter–Few–Rogers covering threshold and use \(B_0=10.79184\) in the 94D fixed-order branch. The count verifies \((2+209/500)B_0^{46}<A_0^{46}\). With \(q=16-B_0\) the fixed determinant is \(p(q-p/4)=530446/28125\) and still exceeds the generic determinant \(3(16-A_0)^2/4\). The predecessor's fixed-order rounding, normal separation and exhaustive two-shear argument apply unchanged.

In 95D use the same fixed-order first two shears with \(B_0\) and put \(D_f=10.55742\). The two unlike-radius merges, finite list of forbidden third-level phases, small nine-simplex replacements and ordinary covering contradiction are all replayed over rationals: \((6+209/500+2/9)D_f^{46}<A_0^{46}\). In the inherited generic branch the determinant stays \(35423551746542630582753447/500000000000000000000000\), smaller than the new fixed determinant \(24961975073/351562500\). The exact density expressions at squared minimum 16 are \(2^{94}\pi^{47}/(47!\,47^7\sqrt{530446/28125})\) and \(2^{143}\pi^{47}/(95!!\,47^7\sqrt{24961975073/351562500})\), respectively. Strict gain uses previous squared determinants \(8488711/450000\) and \(12783236155039/180000000000\), not rounded row comparisons.

Run the adjacent standard-library Node scripts verify-cover.js and verify-tails.js. The radius script is self-contained; the tail script uses the embedded MW92 parent, covering theorem and geometric merging lemmas proved in the linked earlier packages.
