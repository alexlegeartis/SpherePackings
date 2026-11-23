# A 38-dimensional periodic packing from the E7 Hermitian glue

**Exact center density: \(8\sqrt2\approx 11.31370850\)** at squared minimum 8, approximately **19.33%** above the preceding listed general center density 9.4806808308. This construction is periodic and nonlattice, so the lattice-only lower bound is unchanged.

Use the 35-dimensional five-E7 lattice \(L\) with squared minimum 8 and covolume \(2^{31/2}\), documented and proved in the [companion 34–37D construction](../e7_hermitian35_20261123/). Its code alphabet is \(E_7/2E_7\simeq\mathbb F_{64}\oplus\mathbb F_2\), with the quadratic bit \(q(x,e)=q_0(x)+e\), \(q_0(x)=\operatorname{Tr}_{8/2}(x^9)\), and radical \(r=(0,1)\) of least squared norm 6. The code checks are \(\sum x_i=\sum\alpha_i x_i=\sum e_i=0\) for five distinct \(\alpha_i\in\mathbb F_8\). The period of our 38D packing is \(L\oplus 2D_3\), where \(D_3=\{v\in\mathbb Z^3:\sum v_i\text{ even}\}\). The 16 phases are defined below.

## Eight torus shifts with squared clearance at least 6

Choose any missing slope \(\beta\in\mathbb F_8\setminus\{\alpha_i\}\). Let \(B\) be the polar form of \(q_0\) on \(\mathbb F_{64}\) over \(\mathbb F_2\). Choose a set \(S\) of eight vectors, including 0, with \(B(a,b)=1\) whenever \(a,b\ne0\) are distinct. Existence follows by taking the alternating form with Gram \(J+I\) on seven vectors: it has six-dimensional nonsingular quotient by the all-one radical and may be identified with the polar space of \(q_0\). Our certificate checks explicit labels.

For \(a\in S\), choose a shift \(w_a\in E_7^5\) with field syndromes \((s_0,s_1)=(a,\beta a)\) and radical parity \(e=1+q_0(a)\) when \(a\ne0\), or 0 when \(a=0\). Every pairwise difference of distinct shifts has analogous syndromes \((c,\beta c)\), with \(c\ne0\), and its parity is \(1+q_0(c)\). This uses \(B(a,b)=1\) when both labels are nonzero; pairs involving zero hold by construction.

In such a difference class, one active projected slot cannot realize the two syndromes, because \(\beta\) is missing. Three or more active slots have squared length at least 6. With exactly two active slots, the two symbols are \(c t_i,c t_j\) for scalars in \(\mathbb F_8\) summing to 1. Their \(q_0\) values sum to \(q_0(c)\), by the same norm-and-trace identity as in \(L\). They cannot both have norm-2 phases with all other sites zero, which would need radical parity \(q_0(c)\). Thus the squared horizontal distance of distinct shifts modulo \(L\) is at least 6.

## Combining the two parities of the vertical integer lattice

Use the eight labels of \(D_3/2D_3\). In the vertical set \(D_3\), assign its eight residue labels bijectively to the eight shifts \(w_a\). In the other vertical coset \(e_1+D_3\), with \(e_1=(1,0,0)\), use the same eight labels but the shifts \(u+w_a\), with \(u=(r/2,r/2,r/2,r/2,r/2)\). At each vertical point place the entire corresponding horizontal translate of \(L\).

Within one parity, if the vertical labels differ then the vertical squared distance is at least 2, and the horizontal at least 6. Two distinct vertical points of the same label are separated by a nonzero vector of \(2D_3\), of squared length at least 8. In a single fiber the horizontal lattice \(L\) supplies the minimum. For opposite parity, the horizontal difference is in \(u+E_7^5\); doubling it gives a radical class in each of the five blocks, hence squared horizontal distance at least \(5\cdot6/4=15/2\). Vertical squared distance is at least 1 for opposite parities of \(\mathbb Z^3\). These cases prove squared minimum at least 8, attained inside any horizontal fiber.

There are 16 translates of \(L\oplus2D_3\). Its covolume is \(2^{31/2}\cdot16\) and the sphere radius is \(\sqrt2\). The center density is therefore \(16\cdot2^{19}/(16\cdot2^{31/2})=8\sqrt2\).

## Replay

Run python verify35.py to replay the horizontal local alphabet and glue, then python verify38.py to check the explicit eight shifts, all 280 projected two-slot cases for their 28 differences, the vertical basis, covolume, and motif. The exact 16 shifts are in motif.json. Their horizontal coefficient rows are doubled and use the same E7 Cartan basis as certificate.json; their three vertical coordinates are ordinary integers.

The eight nonzero-or-zero field labels have binary span of rank six, so an eight-element list cannot be an additive subgroup. If the packing were a lattice, translating so the origin is a center would make the assignment on the subgroup (D_3/2D_3) a homomorphism to the horizontal torus; projecting its phases to their first field syndrome would then give an additive eight-element image. The independent script `python -B independent38.py` checks this rank and all 28 pair distances via an exact min-plus calculation in a different field model. The comparison is to the [preceding general catalog value](../../LOWER_BOUNDS_1_96.md) from [the 2026 antipode construction](https://arxiv.org/abs/2607.20359).
