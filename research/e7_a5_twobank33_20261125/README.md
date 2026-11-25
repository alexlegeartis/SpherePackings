# 33 dimensions: two A5 banks over four E7 factors

This exact, nonlattice periodic sphere packing has squared minimum distance 8, center density **4/√3**, and ordinary density **4.30349290013662 × 10⁻⁶**. This improves the historical Elkies Q33 center density 3^(33/2)/2^25 by **3.9318248344%**; it does not change the lattice-only bound or establish optimality. See the [Nebe–Sloane historical catalogue, row 33](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html) (Elki95, Elki96a) and [Cohn's historical survey, Table 1](https://dspace.mit.edu/server/api/core/bitstreams/6eda1660-ef8c-4d00-80cf-50805b5b2324/content). The configuration data reuse the E7 Gram matrix and phase isometry from the repository's 35D E7 construction; the 33D product gluing and the two-coset motif are given here.

## Definition

The ambient 33-space is the orthogonal sum of the sum-zero hyperplane in Euclidean R^6 and E7^4 tensor R. Use the seven simple-root coefficient coordinates and Gram matrix in cert.json in each E7 copy. For z in E7^4 compute S(z) blockwise from its coefficient parity masks. Write each mask uniquely as M[t+64e] for t in F64, e in F2, where M is the binary-linear permutation phase_isometry. The block with multiplier a contributes the syndrome (t,at,e), packed as t + 64(at) + 4096e; add syndromes by XOR. Multipliers of the four blocks, in order, are 1,14,15,22. F64 uses modulus X^6+X+1 (integer 67). Write H=ker S and H_s=S⁻¹(s). The kernel H is a full 28-dimensional lattice. S is onto: two distinct multipliers solve the two field equations, and any remaining bit can be set by the radical.

For an integer six-vector z put b(z)=XOR of b_i over its odd coordinates, with b=[749,3755,0,88,403,379]. Set u=3318. Let r be the simple-root coefficient vector with binary mask 82, and T=(0,r/2,r/2,r/2). Take the union of:

1. (z,v) with z integral of sum zero, v in H_{b(z)}.
2. (y−(1/2,...,1/2), T+v) with y integral of sum three, v in H_{u XOR b(y)}.

The first bank is a lattice K. The second is a single translate t+K. They are disjoint by their first five-dimensional coordinates. The doubled offset has local parity mask 63 but exterior syndrome S(2T)=4096, not b(all ones)=3318. Thus 2t is not in K, and the union (which contains zero) is not a lattice.

## Minimum distance

Within either bank differences belong to K. From the first bank to the second they belong to t+K; reversing the order negates the vectors and preserves lengths. In the first case the parity mask has even weight w. The least squared local norm is w (choose equally many +1 and −1). In the second case odd weight w=3 attains 3/2 with three positive and three negative halves; for w=1 or 5 the nearest vector must replace at least one signed half by a signed three-halves to enforce sum zero, which is sufficient, so the minimum is 7/2.

| Weight | Number of masks | Minimum local squared norm | Minimum exterior squared norm |
| --- | ---: | ---: | ---: |
| 0 | 1 | 0 | 0 |
| 1 | 6 | 7/2 | 9/2 |
| 2 | 15 | 2 | 6 |
| 3 | 20 | 3/2 | 13/2 |
| 4 | 15 | 4 | 4 |
| 5 | 6 | 7/2 | 9/2 |
| 6 | 1 | 6 | 6 |

For the zero mask separately, a nonzero local vector in 2A5 has squared norm at least 8. A nonzero vector of H also has squared norm at least 8: all nonzero syndrome-zero phase assignments pass the exhaustive check, and within the all-zero phase the periods 2E7 have minimum 8. A doubled simple root realizes 8.

## Verification, normalization, density

Run `python verify.py` (Python 3 and NumPy). The positive diagonal decomposition of 60G is verified in integers. Recursive enumeration uses integer square roots to enumerate every integral E7 vector of squared norm at most 24, enough to cover every required doubled-coordinate residue mod 4. The checker first recovers the minima for all 128 residues of each block type, and then independently checks the exterior norms for every one of the 64 six-coordinate parity masks. For each syndrome, it enumerates all 2^15 solutions for the four blocks by solving the two field equations and one parity equation. Displayed values in its output are FOUR TIMES the squared norms. There is no floating-point geometric tolerance.

Since det G=2, E7^4 has covolume 4, and H has covolume 4 × 8192 = 32768. The projection K→A5 is surjective with kernel H, hence covol(K)=32768√6. The motif has two points per period, and the packing radius is √2. The center density is therefore 2(√2)^33/(32768√6)=4/√3. Multiplying by π^(33/2)/Γ(35/2) gives ordinary density.

Historical priority is relative to the quoted comparators; no claim of an exhaustive literature search or of third-party endorsement is made.

The independent `node verify_dp.js` replay combines all 8,192 exterior syndromes by integer min-plus convolutions instead of solving individual field systems. It verifies the same nonzero kernel minimum and all 64 local parity classes.
