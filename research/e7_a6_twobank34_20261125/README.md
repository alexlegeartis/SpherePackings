# 34-dimensional two-bank A6/E7 periodic packing

This nonlattice packing has squared minimum distance **8**, center density **8/√7**, and ordinary density **0.00000240447512000689**. It improves on the repository's [34D primitive E7 section](../e7_hermitian35_20261123/) (center 4/√3) by **30.93073414%** and on historical Elkies Q34 (center 3^(33/2)/2^25) by **36.07870126%**. The [historical Nebe–Sloane catalogue, row 34](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html) attributes Q34 to Elkies. This result does not improve the lattice-only row, prove optimality, or claim independent literature exhaustiveness.

## Construction

Use the horizontal lattice H = ker(S) in E7^4 from [the earlier 33D A5 construction](../e7_a5_twobank33_20261125/). All data are supplied independently here. In simple-root coordinates with Gram matrix G from cert.json, represent each E7 coefficient parity mask as M[t+64h] using the binary-linear permutation M in the certificate, with t in F64 and h in F2. The four blocks (in order) have multipliers a=1,14,15,22. A block contributes the binary syndrome (t,at,h), encoded by t+64(at)+4096h, and S is their XOR. F64 uses X^6+X+1. H_s denotes S^{-1}(s). The two distinct multipliers solve the field equations and their radical bits achieve either parity, so S is onto.

Take B=[749,3755,0,88,403,379,5130], u=3318, and let b(y) be the XOR of B_i over the odd coordinates of the integer seven-vector y. For the E7 coefficient vector r with binary mask 82 put T=(0,r/2,r/2,r/2). In the sum-zero hyperplane of R^7 times E7^4 tensor R, take:

1. K: (z,v), with z integral of sum zero and v in H_{b(z)}.
2. t+K: (y−(3/7)1,T+v), with y integral of sum three and v in H_{u XOR b(y)}.

The first bank K is a lattice and the second is one translate. They are disjoint. Their union is not a lattice: the local component of 2t has fractional part −6/7 in every coordinate, which has neither of the two allowed fractional parts 0 and −3/7 modulo one.

## Exact distance proof

Within a bank differences are K; from first to second they are t+K; the opposite differences negate these and have the same norms. For even weight w of the parity mask of z the local minimum is w (choose equal numbers of +1 and −1). For odd weight w of y the local squared norm is |y|²−9/7, and its least possible value is w+2 max(0,3−w)−9/7. For w≥3 use signed units on all odd positions to give sum three. For w=1 one additional even entry of size two is necessary and sufficient (squared norm five). The exterior minimum depends on the parity syndrome, and all masks of each weight give the same value:

| Weight | Masks | Least local squared norm | Least exterior squared norm |
|---:|---:|---:|---:|
| 0 | 1 | 0 | 0 |
| 1 | 7 | 26/7 | 9/2 |
| 2 | 21 | 2 | 6 |
| 3 | 35 | 12/7 | 13/2 |
| 4 | 35 | 4 | 4 |
| 5 | 21 | 26/7 | 9/2 |
| 6 | 7 | 6 | 6 |
| 7 | 1 | 40/7 | 13/2 |

For the zero mask a nonzero local 2A6 vector has squared norm at least eight. In H a nonzero syndrome-zero phase assignment has squared norm at least eight, and the all-zero phase has periods in (2E7)^4, also of squared norm at least eight. A doubled root realizes eight. The least mixed squared distance is 115/14, strictly larger than eight.

Run `python verify.py` (NumPy) and `node verify_dp.js`. Both recheck an exact positive diagonalization of 60G and enumerate integral doubled E7 coordinates with a finite exact squared norm cutoff of 24, sufficient and checked for every needed residue modulo four. The first independently solves the binary field equations and enumerates all 2^15 phases per syndrome; the second instead combines all 8192 syndromes by four min-plus convolutions. No floating-point tolerance enters the minimum certificate.

Finally det G=2, hence covol(H)=4·8192=32768, and K projects onto A6 with kernel H. Thus covol(K)=32768√7. Two centers per cell at packing radius √2 give center density 2(√2)^34/(32768√7)=8/√7; multiply by π^17/17! for the occupied density.
