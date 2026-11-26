# A 14-dimensional periodic packing with density 25/24 times the laminated benchmark

**Center density:** \(25/(384\sqrt3)=0.03758790815036626\ldots\).  
**Occupied-volume density:** \(25\pi^7/(384\sqrt3\,7!)=0.02252510008588239\ldots\).  

This is **4.1667% above** the \(\Lambda_{14}\) center density \(1/(16\sqrt3)\) documented by the [Nebe–Sloane catalogue](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html), and the 14D general bound pinned in this repository from [Cohn's table](https://cohn.mit.edu/sphere-packing/). It is an explicit periodic **nonlattice** construction, not an optimality claim. No lattice-only bound changes.

## Two binary codes on a triangular grid

Take horizontal period \(2\mathbb Z^{12}\). The vertical grid has unit basis vectors e,f with inner product 1/2, so its squared norm at (m,n) is m²+mn+n². Give site (m,n) color c=m−n mod 3. All arithmetic on binary words and five-bit syndromes below is XOR.

- At c=0, place the 144-word binary code C plus the all-half vector, with an additional offset 0, e₀ or e₁ according to n mod 3. Construct C from the [Leech–Sloane description, §2.6](https://neilsloane.com/doc/Me14.pdf): take the all-zero row and 11 cyclic shifts of the residue row supported on {0,1,3,4,5,9} modulo 11 with a prefixed zero. All binary sums of two distinct rows and their complements give the 132 Steiner hexads. Add six disjoint weight-two words {0,1},{2,3},…, {10,11} and their complements. C has minimum Hamming distance 4.
- At c=1,2, place a syndrome coset of the shortened extended Hamming code with the twelve parity-check columns 16+i for 0≤i<12 (five-bit vectors). Each coset has 128 words and minimum distance 4. Let α=28, β=29, γ=30 denote three omitted columns. A coset at (m,n) has label `(m mod 2)(β⊕γ) ⊕ (n mod 2)(α⊕γ) ⊕ [c=2](α⊕β⊕γ)`.

## Separation certificate

In one sheet the binary minimum gives squared distance ≥4, and self-translates by the horizontal period also give minimum 4. Across two sheets with different integer/half-integer types, every horizontal coordinate differs by a half-integer, so the horizontal squared gap is ≥3. The shortest nonzero vertical gap has squared length 1.

For two binary cosets at vertical squared distance 1, the labels differ by α, β or γ: these are the label changes in the three directions (1,0), (0,−1), (−1,1) from color 1 to color 2. Such a syndrome has no word of weight one (it is an omitted column), or of weight two (it has odd parity in the fifth bit). The horizontal squared gap is therefore at least 3.

At vertical squared distance 3, the third color difference is 0. For binary sheets the label differences are three distinct nonzero even syndromes β⊕γ, α⊕γ, α⊕β. For half-integer sheets the offset changes, because the six vertical directions of squared norm 3 change n nontrivially modulo 3. The differences between two offsets have Hamming weight 1 or 2, so two translated C sets are disjoint by its distance 4. In either case distinct horizontal parity cosets give squared horizontal gap at least 1. Larger vertical squared distances are already ≥4. The true minimum is exactly 4 by the horizontal period, so the spheres can have radius 1.

The periods 6e and 6f preserve all labels and offsets. There are 12 sheets of each color in this cell and \(12(144+128+128)=4800\) sites per horizontal/vertical product cell. Its volume is \(2^{12}\cdot36\sqrt3/2=73728\sqrt3\). The ratio gives the density above. This set cannot be a lattice or a translate of one: fibers in the vertical projection would then be cosets of a common horizontal lattice and have the same number of sites per horizontal period, whereas these fibers have 144 or 128.

Run `python3 verify.py` for a standalone finite integer replay of the 36 sheets, all short vertical offsets and every short horizontal pairing. Source attribution: the 144-word code is historical; the new ingredient is the mixed triangular stack, in which the three honeycomb edge directions use three distinct missing Hamming columns.
