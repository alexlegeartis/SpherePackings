# Local phase obstruction for two published MW80 sections

This is an exact local result for fixed periods. It neither establishes a
global optimum nor improves any numerical packing bound. In dimension 79 it
applies to the two-translate MW80 construction in the live density table.
In dimension 78 it concerns the three-translate MW80 construction in
creative_sections, not the stronger live Γ72-based entry.

## Statement

For the 79-dimensional construction, every sufficiently small nonzero change
of the relative translation strictly decreases the minimum squared distance,
if the section lattice is held fixed. In the 78-dimensional construction the
minimum cannot increase in a neighborhood, even when all three translations
are permitted to move independently. If the relative translation of the
limiting pair changes nontrivially, its minimum strictly decreases. A common
translation is a symmetry.

## Proof

Use the row-coordinate convention and the parent Gram matrix G from the
[existing exact antipode construction](https://github.com/alexlegeartis/SpherePackings/blob/6696a479a8a55e440ff8e65f10702702fe7d599c/research/creative_sections_20261116/PROOF.md).
The underlying Mordell–Weil parent and minimum 14 are as in Schütt–Shioda,
*Mordell–Weil Lattices*, Theorem 13.35, and in the linked exact reconstruction.
The input normal and cluster labels are the same as that construction.

Write B for its normal-label matrix, K=B G^(-1) B^T, and d for the
difference of labels in the limiting pair. The orthogonal projection of a
parent vector x with label d is v=x−d K^(-1) B G^(-1). The published
projection-loss proof gives the minimum squared distance m=14−β, with the
attainment vector at that pair. The verifier recomputes β over all label
differences: it is 779/2572 in rank 78 and 41/234 in rank 79.
Each literal x in the new certificate has label d and squared norm exactly
14. Its projected vector v is therefore an actual minimum contact from
one translate to the other, not a putative nearest-image direction.

In each dimension the supplied projected contact rows span the entire
horizontal space; rank modulo 101 after clearing denominators certifies
this directly. There are strictly positive rational coefficients a_i with
sum one and sum a_i v_i=0. The verifier reconstructs the coefficients from
a small prescribed support plus a common uniform positive coefficient by
exact linear equations. It checks all coordinates and positivity. Hence
for each unit horizontal direction u one of these contacts satisfies
v_i·u<0: if all were nonnegative then positivity and the balance would
make all zero, contradicting span. Compactness of the unit sphere gives
a uniform ε>0 with min_i v_i·u≤−ε. If the relative translation of the
limiting pair changes by w, the same literal contact differences move
to v_i+w. Therefore their minimum squared norm is at most

    m − 2ε ||w|| + ||w||².

It is strictly lower than m for sufficiently small nonzero w. When the
limiting relative translation does not change, its contacts remain,
so the overall minimum is still bounded above by m. This proves both
claims, with no need to enumerate all additional short vectors or assume
that the other translation differences stay close to the original normal
projection labels. It makes no claim about perturbing the period Gram,
about more translates, or about distant choices of the same phases.

## Verification and scope

Run OPENBLAS_NUM_THREADS=1 python -B verify.py (NumPy and SymPy).
The script uses the pinned source parent, its exact integer inverse after
multiplication by 41, integer Gram products and exact rational equations.
The floating inverse is only a proposal checked by an exact multiplication.
For period primitiveness, the parent minimum and the original construction
see the prior proof above. The included source parent is from
[cm58 frontier/parent80.npz](https://github.com/alexlegeartis/SpherePackings/blob/6696a479a8a55e440ff8e65f10702702fe7d599c/research/cm58_20261116/frontier/parent80.npz)
(Git blob f62fa88c44f4ce90563e6e5d0219a7bf5769f30d). The same inputs can be
obtained from those pinned repository paths. The checked set of literal
contacts suffices; no assumption of completeness or reliance on the
search's LP status enters the proof.
