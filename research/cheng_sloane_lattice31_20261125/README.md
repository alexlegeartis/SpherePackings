# A 31-dimensional lattice of center density √2

This explicit additive lattice attains center density √2 and occupied density at least 1.38411060495194e−5. It improves the preceding [catalogued lattice bound Q31](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html), of center density 3^15/2^(47/2), by a factor 2^24/3^15 ≈ 1.1692330. It **ties, but does not improve**, the repository's [31D unrestricted lower bound](../vardy30_shear31_20261123/). This is a reproducible construction, not a classification or a claim of optimality or exhaustive historical priority.

## Construction

Take the binary [32,17,8] code of [Y. Cheng and N. J. A. Sloane (1989)](https://doi.org/10.1137/0402003), with the explicit generator as displayed in [M. Grassl's code tables](https://www.codetables.de/BKLC/BKLC.php?q=2&n=32&k=17). Shorten on coordinate zero: retain words zero there, remove that coordinate and index the remaining coordinates from zero. The resulting code C is even, of dimension 16 and minimum distance 8. Both binary generator matrices are pinned as little-endian integer masks in witness.json. The authors also describe a 32D packing of center density 2; that is distinct from the claim here.

Let S={0,...,6}, T={7,...,30}, and let t be the 0/1 indicator of T. Exact binary enumeration verifies:
* t belongs to C;
* each projection of C to S has even weight;
* every c in C supported on T has weight divisible by four.

Let D31 be the integer vectors of even coordinate sum and e₀ the first coordinate vector. Set q(c)=wt(c)/2 mod 2. Define
```
K = union over c in C of (2c + 4q(c)e₀ + 4D31),
z = t + 4e₀,                    L = K ∪ (z + K).
```
Here binary words in real vectors are represented by 0/1 coordinates, not by signed values.

**Additivity.** q(c⊕d)=q(c)+q(d)+c·d mod 2. Passing from ordinary c+d to c⊕d in 2c+2d costs 4(c∧d), whose carry parity is c·d. Thus K is a lattice, with 2^16 distinct classes modulo 4D31. Since q(t)=0, 2z=2t+8e₀ lies in K, whereas z is not even, so L is a lattice with 2^17 classes.

**Minimum.** In a nontrivial binary residue of K at least eight coordinates are 2 modulo 4, giving squared norm at least 32. The pure period 4D31 has squared minimum 32. For z+2c+4q(c)e₀+4D31 the 24 coordinates on T are odd. If c is nonzero on S, at least two of those exterior coordinates are 2 modulo 4, adding at least eight. If c vanishes on S, then q(c)=0. To realize the nominal minimum 24, all exterior coordinates must be zero and every T-coordinate must be ±1. Reducing t+2c to these nearest values subtracts 4 at exactly wt(c) positions, an even number. The extra 4e₀ enforces the opposite total carry parity. Hence the nearest choices cannot be taken together: an odd coordinate must have absolute value ≥3, or an exterior coordinate absolute value ≥4. This raises the squared norm by at least eight. The period contains 4(e₀+e₁), so 32 is attained.

**Density.** The covolume is (2·4^31)/2^17=2^46, Gram determinant 2^92, and packing radius √8. The center density squared is 8^31/2^92=2. Occupied density is this center density times the 31D unit-ball volume. The unrestricted construction has exactly the same center density; no unrestricted record changes.

## Replay

Run `python3 verify.py` and `node independent.js`. The first checks the binary source, exact hypotheses, index and comparison. The independent verifier generates the lattice directly from integer vectors; it verifies 131,072 different residue classes, their closure under generator doubling, and the exact minimum in every nonzero class by a two-state distance dynamic program for the parity of all integer carries. Data and computed summary: witness.json and results.json. No downloaded dataset or numerical tolerance is needed.
