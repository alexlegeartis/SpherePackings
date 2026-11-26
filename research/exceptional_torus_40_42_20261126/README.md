# Wrapped antipode anticodes on the exceptional normal sections in 40–42D

No record. For ANY even unimodular 48D lattice of minimum 6 admitting the indicated primitive normal K=√3 E_r (r=8,7,6 respectively), the centre density certified by the standard normal-diameter antipode formula on this K never improves its bare one-class section. This is for all cardinalities, all radii and all wrapping of normal representatives. Unlike the 49D theorem this result pertains to the normal-diameter certificate, not to the true horizontal distances of arbitrary discriminant motifs. We actually assembled a nine-class candidate of exact squared torus diameter 2/3 in each dimension; nine is sharp at that diameter, but improving would require 11 in 40D and 12 in 41–42D. All larger diameters were considered, not merely that threshold. The full singleton densities are δ40=(3/2)^20/81, δ41=(3/2)^(41/2)/sqrt(2·3^7), δ42=(3/2)^21/sqrt(3^7), the corresponding public floors at https://github.com/alexlegeartis/SpherePackings/blob/98b73324217feb9914b0dd1eea4c25aacf90ed93/LOWER_BOUNDS_1_96.md . Classical construction background: Conway–Sloane, Laminated lattices, Corollary 8; section notation follows https://arxiv.org/abs/2607.20359 .

# Reduction to an anticode

Let C be the simple-root Gram, d=det C, m=3d. Under scaling by √3 the group of canonical classes is E_r^*/3E_r (size d·3^r). For integer weight labels c in fundamental-weight coordinates, class key is adj(C)c mod m and torus squared norm is the minimum of n/m where n=cᵀadj(C)c over representatives. If S has M points and squared torus diameter β=b/m, density ratio to the singleton is M(1−β/6)^((48−r)/2). A motif meeting a target radius always admits at least its diameter loss in this sense; a further radius reduction cannot help. The entire problem here is whether this ratio may exceed one. The structural move is to find a **long-difference code** H in the same finite group. When every nonzero word of H has square norm >β, S uses at most one point of every H-coset. Thus M≤d·3^r/|H|, independently of winding, without any semidefinite or numerical optimization.

# Small diameters

For β≤2/3 choose shortest lifts from a common packet point, taken as origin. An edge which does not lift coherently gives a nonzero K-period of length ≤3√β≤√6. Equality throughout forces the three collinear equally spaced representatives 0,a,−a with 3a a shortest K-period. For any different fourth point z, its edges to a and −a cannot wrap, since equality in the same triangle inequality would force z=−a or a. Distances ≤|a| to both imply z=0 by the parallelogram law. Hence packets of at least four sites lift coherently even at equality. Partition these lifts by the root-lattice coset in E_r^*/E_r after rescaling. Within one bank any nonzero difference at β≤2/3 is a root of squared length 2 in unscaled units. Cross-bank distances (within the ≤2 threshold) are all equal: 4/3 for E6 or 3/2 for E7. This follows also from the small norm table and the fact that every norm-2 weight is in the root lattice (nonzero discriminant classes have nonintegral even norm residue; verify with the Cartan inverse). Each bank is a regular simplex; the affine direction spaces of two banks are orthogonal by subtracting four equal cross-distance squares. Therefore M≤r+#banks≤9. Below 2/3 there is only one site per bank, giving E6 cap3, E7 cap2. E8 itself has one bank. The stated nine-site constructions at threshold are given in the replay.

# Integer upper gates

Each line is b -> maximum M used; the last general cap applies to all larger diameters. For E7/E8 integer b=13/7 respectively need not actually occur, but the safe general cap already wins there. Each cap below gives STRICTLY smaller density than the singleton, checked by M²(6m−b)^(48−r)<(6m)^(48−r).

- E6: 4→3, 6→9, 10→27, 12→81, 16→81, ≥18→2187.
- E7: 3→2, 4→9, 7→54, 8→162, 11→162, 12→2187, ≥13→4374.
- E8: 2→9, 4→81, 6→2187, ≥7→6561.

The large-diameter caps (other than full-group and E7 final involution) come from six supplied ternary row spans in root coordinates; each is checked to avoid every forbidden difference. Compute all distances through β=2 by integer enumeration of coefficients |c_i|≤√12<4, by testing c against the ith simple root. All determinant and norm arithmetic within Number in the small dimensions is integer and exact; exact density powers use BigInt. No float-based SVP, LP or graph optimizer is required. Defining the Cartan graph as a chain 0…r−2 plus an edge from 2 to r−1 fixes all coordinate conventions.


Run `node verify.js`. The script checks the six exact separating codes and three constructive nine-class motifs. The small-diameter structural bound is analytic above.
