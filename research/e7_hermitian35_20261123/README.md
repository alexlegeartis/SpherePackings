# Four sphere-packing constructions from five-block E7 Hermitian glue

All four lattices have exact squared minimum 8.

| Dimension | Gram determinant | Center density | Preceding listed center density |
|---:|---:|---:|---:|
| 34 | \(3\cdot2^{30}\) | \(4/\sqrt3\) | 2.2220345697 |
| 35 | \(2^{31}\) | \(4\) | \(2\sqrt2\) |
| 36 | \(2^{30}\) | \(8\) | 5.2816290013 |
| 37 | \(2^{31}\) | \(8\) | \(4\sqrt2\) |

Center density omits the unit-ball volume. These are lattice and hence general periodic lower bounds, not optimality claims. The comparison rows are the preceding catalogue values.

## Local alphabet

Take the standard E7 Cartan Gram, diagonal 2 and off-diagonal −1 on the tree with arm lengths 1, 2, and 3. It is positive definite with determinant 2 and squared minimum 2. On \(E_7/2E_7\) use the quadratic form \(q(v)=\|v\|^2/2\bmod2\). Its radical for the polar form has order 2. Write \(r\) for the nonzero radical class or its representative with coefficients 1 on the three nonadjacent vertices (mask 82) and 0 elsewhere. One has \(q(r)=1\) and least squared norm 6: the representative has squared norm 6, and exact enumeration of all vectors of squared norms 2 and 4 shows none in that class. A nonzero phase of \(q=0\) has least squared norm at least 4, and every nonzero phase costs at least 2.

Identify the binary phase space with \(\mathbb F_{64}\oplus\mathbb F_2\) so \(r=(0,1)\) and
\[
q(x,e)=q_0(x)+e,\qquad
q_0(x)=\operatorname{Tr}_{\mathbb F_8/\mathbb F_2}(x^9).
\]
Here \(x^9=\operatorname{Norm}_{\mathbb F_{64}/\mathbb F_8}(x)\). The polar form of \(q_0\) is nondegenerate: fixing \(x\ne0\), its polar value on \(y=xz\) is \(\operatorname{Tr}(x^9(z+z^8))\). The map \(z\mapsto z+z^8\) surjects onto \(\mathbb F_8\) and multiplication by \(x^9\) is invertible there. Fix a polar isometry of the quotient, then shift each lifted basis vector by \(r\) if necessary to match the quadratic terms, since adding \(r\) toggles \(q\). The exact 128-entry isometry appears in certificate.json and is replayed by verify.py.

## The 35D lattice

Choose five distinct \(\alpha_i\in\mathbb F_8\), and define the phase code by
\[
\sum_i x_i=\sum_i \alpha_i x_i=0,\qquad \sum_i e_i=0 .
\]
The two field checks and one parity check give binary rank \(3\cdot6+4=22\). Let \(L_{35}\subset E_7^5\) be the inverse image of this code modulo \(2E_7^5\).

A nonzero projected \(x\) has at least three occupied slots. At exactly three slots, with the chosen \(\alpha_i\) distinct, its symbols have the form \(z(c_i,c_j,c_k)\), where \(c_i,c_j,c_k\in\mathbb F_8\) and \(c_i+c_j+c_k=0\). Hence
\[
\sum_\ell q_0(x_\ell)=\operatorname{Tr}\left(z^9\sum_\ell c_\ell^2\right)=0.
\]
If those three phases all cost 2 and the other slots are zero, the sum of the quadratic bits is 1, incompatible with even \(e\)-parity and the displayed identity. So a weight-three \(x\) costs at least 8: either some active slot costs at least 4, or a dormant slot is a nonzero radical and itself costs at least 6. Four or more active slots already cost at least 8. If \(x=0\), a nonzero phase has at least two nonzero radical slots by parity, costing at least 12. The zero phase consists of the period \(2E_7^5\), whose nonzero vectors have squared minimum 8. Thus the full squared minimum is 8.

The period has covolume \(2^{35}2^{5/2}\) and index \(2^{22}\) in \(L_{35}\), giving lattice covolume \(2^{31/2}\).

## The 34D section

In the first block take the fundamental weight dual to the last end vertex of the longest arm (index 6). It has squared norm \(3/2\), directly from the inverse E7 Gram. It belongs to the dual of \(L_{35}\) because \(L_{35}\subset E_7^5\), and it is primitive because some glue generator pairs oddly with it and the period includes a vector pairing to 2. The perpendicular section has covolume \(2^{31/2}\sqrt{3/2}\) and squared minimum 8, attained by a doubled root in any other block. This yields center density \(4/\sqrt3\).

## An order-four extension to 36D

Let \(u=(r/2,r/2,r/2,r/2,r/2)\) in the five E7 blocks. Stack \(L_{35}+ku\) at heights \(k/\sqrt2\), \(k\in\mathbb Z\), to form \(L_{36}\). As \(4u\in2E_7^5\), the horizontal phase depends only on \(k\bmod4\). For odd \(k\), doubling any vector of \(L_{35}+ku\) gives a radical residue in each local block, so the original horizontal squared norm is at least \(5\cdot6/4=15/2\), and the vertical squared norm is at least \(1/2\). For \(k=2\bmod4\), the horizontal phase at \(x=0\) has odd radical parity and costs at least 6; with \(x\ne0\) at least three active slots cost at least 6. The vertical part adds at least 2. For nonzero \(k=0\bmod4\) the vertical squared norm is at least 8. The zero layer has squared minimum 8, and the covolume is \(2^{31/2}/\sqrt2=2^{15}\), giving center density 8.

## An order-two extension to 37D

Choose \(w\in E_7^5\) with phase syndromes \(\sum x_i=0\), \(\sum\alpha_i x_i\ne0\), and odd total \(e\)-parity. The supplied vector uses \(x_0=x_1=1\), other \(x_i=0\), and only \(e_0=1\). In \(w+L_{35}\) every phase has at least two active projected slots. If exactly two are active their symbols are equal, so their \(q_0\) values coincide. Thus two norm-2 phases with other slots zero would have even \(e\)-parity, a contradiction. The squared depth of this coset is at least 6. The coset \(w+2u+L_{35}\) still has at least two active slots and therefore squared depth at least 4.

The shift \((w,0)+L_{36}\) has squared depth at least 6: the odd old layers have horizontal squared length at least \(15/2\) by the half-radical argument, the \(k=2\bmod4\) old layers combine horizontal at least 4 and vertical at least 2, the zero old layer has at least 6, and nonzero old \(k=0\bmod4\) adds vertical at least 8. Also \((2w,0)\in L_{36}\). Adjoin a new layer generator \(((w,0),\sqrt2)\). Odd new layers then have squared length at least \(6+2\), even nonzero new layers have at least 8 from the newest coordinate alone, and the zero layer inherits the old minimum. Its covolume is \(2^{15}\sqrt2\), yielding density 8.

## Exact replay

Run python verify.py and python basis.py (Python 3 and SymPy). The former checks the exact local norm-2/4 census by inverse-Gram coordinate bounds, local quadratic isometry, phase-code rank, and every one of the \(64^3\) projected words. The latter constructs the explicit Hermite basis and all four Gram determinants exactly. bases.json stores each basis in coordinates relative to five E7 simple-root blocks followed by the indicated vertical axes. The coefficient rows are multiplied by 2; copies of E7_gram in certificate.json set the horizontal metric and the final listed independent diagonal entries give the vertical metric. The optional --write flag regenerates the corresponding JSON.

## Exact scope and comparison

For fixed period `(2 E7)^5`, the 35D construction is optimal among additive subgroups of the binary phase quotient. At squared minimum 8, an exact nonnegative-Fourier test function bounds even arbitrary nonlinear motifs by 6,936,440 phases. An additive motif has size a power of two, so it has at most 2^22 phases, attained here. At squared minimum 6, a 64-point anticode in one block gives at most 2^29 phases, which has smaller density. At minimum 4 or below, the total number 2^35 of available labels already gives smaller density. This claim concerns this fixed period and binary shifts only.

The independent scripts `e7_direct.py`, `phase_syndrome.py` and `verify_additive_optimum.py` use standard-library integer or rational arithmetic. They check the local alphabet in different coordinates, the exact 37D phase depths in a different field representation, and the 35D Fourier certificate. `results.json` records the exact determinant and squared-density comparisons against the repository's [pinned general](../lower_bounds_1_96_20261122/cohn_inputs.json) and [lattice-only](../lattice_lower_bounds_1_48_20261122/cohn_inputs.json) baselines, including the preceding [36D metric deformation](../../constructions/dim36/).
