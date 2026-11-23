# Exact global exclusion for arbitrary integral four-Leech motifs retaining Γ72

No new lower bound. This is an upper bound for one structured class, not a sphere-packing upper bound in dimension 96.

## Scope

Let Λ be the Leech lattice with raw squared minimum 4. For the first three Λ blocks use the physical Gram factor 1/2 and retain Γ72 as a period lattice. Its glue is characterized by common A-parts and total B-part zero, via the standard totally singular complementary polarization Λ/2Λ=A+B. Use a fourth Λ block with an arbitrary positive squared multiplier (s) and with raw period (2Λ). The packing may take **any subset** of residues
[
  (Λ^3/Γ_{72})\times(Λ/2Λ).
]
It may take any number of horizontal shifts at a vertical residue and use arbitrary nonlinear labels, B-to-A mixing, and nonconstant horizontal fibers. No additive structure is assumed. Horizontal phases outside the integral Λ³, other vertical period lattices, or anisotropic vertical deformation are not covered.

The center density of such a packing is at most
[
  2^{48}(4/3)^{12}
  =4722366482869645213696/531441
  \quad (<0.836\ \text{times the pinned old 96D periodic comparator}).
]
Center density comparisons are made via exact squared rationals.

## Exact obstruction at squared minimum (μ>6)

A horizontal coset has syndrome ((a,b,t)), with (a,b) the first and second A-differences from the third block, and (t) the total B-label. Its torus squared depth is integral and at most six. If (a,b) are dependent, realize the syndrome by a single raw class of squared norm at most eight. Otherwise realize it in two blocks by choosing a raw norm-four class for the first and any class of raw norm at most eight for the second. The Leech parity facts used here are: all nonzero A labels have 24 q4 classes; every anisotropic class has q6; all remaining nonzero classes have q8. The pinned upstream public Γ72 verification (dim73-82-gamma72-depth6-lattice) supplies the complete exact norm-four table and proves these parity facts from the unimodular Leech theta series and the parallelogram law.

Consequently there is at most one horizontal offset per selected vertical residue if (μ>6). In the case (μ-4s>5), every pair of selected vertical residues with raw distance four must have horizontal distance at least six. This implies that their differing horizontal (a,b) define a plane with three nonzero distinct labels. For each pair of those labels (x,y), include a norm-four first occupied block with A label (x) and B label (u) from its 24-element fibre (S_x). The second occupied block then has label ((y,t+u)). It must not represent norm six, since a raw (4+6) vector has physical squared norm five. Thus (y\cdot u=y\cdot t) must be constant over (S_x), consistent when the two possible (x)'s are switched, and linear in (y) across the plane. The verifier enumerates precisely this necessary condition and obtains 975 planes, i.e. 5,850 ordered pairs. It does not require sufficiency or each plane admitting such fifth-coordinate assignments.

The code constructs the Cayley graph of q4 vertical residue differences on (2^{24}) vertices. Exact Walsh transform of its indicator gives degree 98,280 and least eigenvalue (-24). The corresponding transform of the 5,850 ordered pairs has minimum (-390) on 20,280 dual characters. The normalized real vector of these characters at a horizontal ordered A-label ((a,b)) has correlation exactly (-1352/20280=-1/15) across **each** allowed difference, by a second exact Walsh transform.

Write (N) for the occupied vertical residues and (E) for their number of ordered q4 edges. Since the constant eigenvalue is 98,280 and all other eigenvalues at least (-24), PSD gives (E\ge(3/512)N^2-24N). Apply the same lower spectrum bound to the Gram matrix of the unit character vectors (set each omitted residue's vector to zero). This gives (-E/15\ge-24N), and hence (N\le65536).

The center-density expression is (N(μ/4)^{48}/(2^{24}s^{12})). Keeping full periods implies (μ\le\min(8,16s)), so in this branch the upper expression is maximized at (s=1/2) and is at most (2^{52}). In the other branch (μ-4s\le5), take (N\le2^{24}) and (μ\le\min(8,16s,5+4s)). The upper expression increases on (0<s\le5/12), still increases on (5/12\le s\le3/4) (log derivative (12(12s-5)/(s(5+4s))\)), and decreases thereafter. Its maximum is (2^{48}(4/3)^{12}).

## Low actual minima and arbitrary fiber multiplicity

For any vertical residue fiber of horizontal offsets, with squared minimum (μ\le6):
- If (μ\le3), the total (2^{36}) horizontal cosets give the direct bound.
- If (3<μ\le4), split into (2^{12}) slices of constant (a). Within each slice any difference ((b,t)) must be singular for the split quadratic form (b\cdot t\), since anisotropy gives a single raw norm-six block. After translating one member to zero, pairwise singular differences imply the linear span is totally singular, hence has rank at most 12. This gives (2^{24}) horizontal cosets per vertical residue.
- If (4<μ\le5), a repeated (a) gives a single-block representative of raw norm at most eight and is impossible. Thus at most (2^{12}).
- If (5<μ\le6), all pairs must have differences among the same 5,850 ordered pairs just analyzed. Unit character vectors with correlation (-1/15) pairwise have cardinality at most 16 by nonnegativity of the squared norm of their sum.

With (2^{24}) available vertical residues, use (s\ge μ/16). In each interval the bound then increases with (μ), and all four endpoint bounds are strictly below the quoted ceiling. This covers multiple fibers without relying on the solution of any general periodic packing problem.

## Strict replay

Run `python -B verify.py` with a C++20 compiler `g++`. The exact input is the parity and pairing export from the upstream pinned Leech shell verifier, copied from its independent exact replay. Its digest and all regenerated finite counts are checked. The two 24-bit Walsh transforms use signed integers, and the final comparisons use rational arithmetic. The old exact periodic denominator comparator is included literally. The exploratory siblings in the parent directory are not required.
