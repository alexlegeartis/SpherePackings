# A repaired Eisenstein shear in dimension 95

A lattice packing exists with squared minimum 16 and occupied density at least **2.115959060080788e-21**. This improves the previous verified MW92 lattice bound **2.104056220143222e-21** by **0.5657092%** (rounded down). It improves both the lattice and unrestricted lower bounds in the repository. The shifts are existential, not published explicit coordinates. The earlier derived external comparison with the minimum-matched MW92 × A3 product is unchanged. No global optimality or claim of an exhaustive priority review is made.

Use the [existing MW92 order-three dichotomy](../mw92_order3_dichotomy_20261122/) and its [sharpened covering constants and prior 95D result](../mw92_cfr_taylor13_20261123/). Exact comparison was made against the [pinned prior 95D row](https://github.com/alexlegeartis/SpherePackings/blob/ee34b698d529020944962f17ff60e67600463ba2/research/mw92_cfr_taylor13_20261123/result.json). The foundation is Shioda's 92-dimensional Mordell–Weil lattice L with squared minimum 16, covolume 47^7, and integral isometry J satisfying J²+J+I=0; see the [parent and arithmetic audits](../mw92_extensions_20261122/README.md#parent-and-covering-input). In particular this proof uses the Coxeter–Few–Rogers equal-ball covering theorem [DOI](https://doi.org/10.1112/S0025579300002059), not a conjectural unequal-radius version.

## Why fewer balls suffice

Write d for torus distance from zero, in length units; all named depths below refer to its square. For any k specified torus centers, the covering theorem supplies a common shift outside balls of squared radius D_k when k D_k^46 < A^46 with A=11.00098526, since the predecessor certified the one-ball covering-volume inequality strictly at A. Set D6=10.58, D7=10.545, D9=10.487. Let a have d(a)^2>A, as in the preceding construction, and set p0=16−A, t=16−3p0>1; write ω=−J, so ω²=ω−1. Multiplication by each of the six Eisenstein units preserves torus distance, and N(x+yω)=x²+xy+y².

If d((1+ω)a)^2<t, use **exactly the preceding fixed-order-three branch B**. Its fixed class of depth at least 34/3 and the two covering shears give the squared tail covolume 652349944595359843/9000000000000000, as independently recomputed by the checker. In the rest of the proof d((1+ω)a)^2≥t. In particular every norm-one label times a has depth >A, and every norm-three label times a has depth ≥t.

Take three vertical basis directions, the first two with tail Gram 5 times the unit triangular Gram. For coefficient vector (i,j,k), the horizontal translation is (i+jω)a + kb modulo L. A layer k=1 therefore has horizontal squared distance d(b+(i+jω)a)^2. We write integer pairs as Eisenstein labels. In addition to the six nearby labels
```
P=(0,0),(1,0),(0,1),(-1,1),(1,-1),(1,1)
```
only six outer labels need special attention:
```
F=(0,-1),(-1,0),(-1,2),(0,2),(2,-1),(2,0).
```
They lie at squared label distances 7/3 from (1/3,1/3); P lies at 1/3 or 4/3. The actual integer tail check below, not a continuum diagram, guarantees that this list suffices.

**Near a small torsion class.** If d(3a)≤2, reduce P∪F modulo the scalar 3 to nine residues. Choose one representative of each, taking all six P labels first. If d(αa)≤2 for either norm-seven direction α=2+ω or 1+2ω, instead use the seven residues modulo (α). Explicit reductions are x+wy modulo 7, where w is the root 3 or 5 of w²−w+1 annihilating α. In each case representatives include P and every other outer label differs from its representative by a unit times 3, respectively α. Indeed the maximum norm of any difference in P∪F is 9, and the residue/difference statements are checked explicitly. Select b to avoid the negatives of the 9 or 7 representative phases at squared depth D9 or D7. All nearby labels are protected to Dk. Every other outer label has distance at least sqrt(Dk)−2>1. This uses only the triangle inequality.

**Away from those small torsion classes.** Now d(3a)>2 and both norm-seven directions have distance >2. Choose b against just the six negatives of P, to squared depth D6. At most one outer label can have distance below 1. To see why, the nonzero differences between outer labels have norms 1,4,7,9. Norm 1 already gives distance above 2, as do 7 and 9 by the present assumptions. If c is any bad outer label, there is some e in P with N(c−e)=4. Therefore d(2a)≥sqrt(D6)−d(b+ca)>2, since every norm-four label is twice a unit. Two bad labels would have difference of distance less than 2, impossible. The six labels of F are transitive under affine permutations of the nearby unit triangle. We may thus place the bad label at (0,−1) for the tail metric. This is just integer relabeling in the proof, not an assumption of complex-conjugation symmetry of L: for m↦Bm+v preserving the triangle, use the unimodular coefficient transformation (m,k)↦(Bm+kv,k). For k=0 the relabeling preserves the norm-one and norm-three bounds, and for k=1 it permutes P and F. Coefficients at other k require no horizontal bound. If there is no bad label, any orientation works.

## Three exact tail metrics

Use the Gram
```
[5,   5/2, -u]
[5/2, 5,   -u]
[-u,  -u,   c]
```
with (u,c)=(5/2,16−Dk) in the k=7 and k=9 cases. In the six-ball case use the tilt (u,c)=(D6/4,11−D6/2), in the selected orientation. For every tail vector of squared norm q<16, add the proven horizontal floor for its label: in layer zero it is A for norm-one labels and t for norm-three labels, zero otherwise; in layer one it is Dk on P, 1 on F apart from the one allowed exception for k=6, and zero otherwise. Reverse all signs when appropriate. Every sum is at least 16. The exact verifier enumerates the complete coefficient box given by z_i² ≤ 16(G⁻¹)_ii and handles all other layers without a horizontal assumption. It also checks the transitivity, the two norm-seven orbits and residue coverage.

Each packing is the lattice generated by (L,0) and the three paired shear generators. The determinant of its Gram is 47^14 times the determinant of the tail Gram, independently of the actual continuous shifts. The squared minimum is exactly 16 because it is attained in L. The repaired generic-branch tail determinants are 576659/8000, 2273/32, and 11539/160; all are below the unchanged fixed-order branch's 72.4833271772622…. Consequently 4^95/(47^14 times this worst determinant) is a certified **lower bound** for squared center density. The prior worst determinant was 73.3057365619806…. Taking the unit-ball volume into account gives the occupied-volume density stated above.

Run `python3 verify.py` (standard library only, assertions enabled). It prints `result.json`, including a floor on occupied density computed with a rational lower interval for π. The inherited torus covering theorem and previous fixed-order class argument are reused, not independently reproved by that computation.
