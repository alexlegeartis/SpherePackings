# Exact 48D optimum for orthogonal odd-Leech half-period products

Let O be the odd Leech lattice, with covolume 1 and squared minimum 3. Consider any finite motif, not necessarily a subgroup, of half-period translates of H=√2 O ⊕ √(2r) O for real r≥1. Its center density is at most (3/2)^24. P48p attains this, from the two odd-Leech eigensublattices; twice any lattice point is the sum of its integral eigenspace doubles. See [Nebe, Theorem 4.10 and Table 2](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/DIM48.pdf). This is a statement about this family, not a general 48D optimum.

## Normalization and density

Write μ≤6 for the squared minimum (the first period block gives this ceiling), N for the motif size, and q for sixteen times the squared depth of a half-period label in the shorter factor. The ratio to P48 is N(μ/6)^24/(2^24 r^12). The exact odd-Leech depth list is q in {0,24,32,40,48,56,64,72,80}. Therefore if μ>5, the second half-period labels are distinct, and N≤2^24. This gives the asserted optimum. If μ≤4.5, fixing the first label leaves at most (8r/μ)^12 points by [24D Leech optimality](https://arxiv.org/abs/1603.06518). The density ratio is then at most (2μ/9)^12≤1.

For 4.5<μ≤5 use the following real functions on O/2O. Each has nonnegative discrete Fourier transform. The first has origin value 2048 and value 1 at q24, zero elsewhere; that graph has valency 2048 so Fourier positivity is immediate. The second has origin 72 and value 1 at q32, zero elsewhere. The third has the respective values
```
13210188, 205176, 86860, 16512, 2448, 7032, 6732, 0, 0.
```
The replay checks positivity of the latter two on all 2^24 characters. For any motif with no nonzero differences in a function's support, its size is at most the ratio of origin value to mean, by Fourier positivity applied to the sum over all ordered motif pairs. These ratios for the three functions are C1=2^23, C2=72·2^24/49200, C3=4403396/1559.

* If r≥27/20, fix the second label. Within each fiber the first-label minimum exceeds 4.5, so the third function gives N≤2^24 C3. This is below the count needed even at μ=5 and r=27/20.
* If smaller r and μ>2+2r, the product of the second function in both factors vanishes on every nonzero motif difference. Thus N≤C2², below the count needed even at μ=5 and r=1.
* Otherwise use the second function in the first factor and the first in the second. Its support has squared distance at most 2+1.5r<4.5. The resulting cap C2 C1 is below the required N≥2^24 r^12 (6/(2+2r))^24, which is decreasing for r≥1 and hence smallest at r=27/20.

All comparisons are strict in this middle range and checked rationally. In the low range equality μ=4.5 would require an integral fiber to have nonintegral (8r/μ)^12 points at r=1. In the high range equality requires r=1, μ=6 and N=2^24.

## Exact alphabet and Fourier replay

In physical coordinates scaled by √8, construct O from the even part {2c+4h: c in the extended binary Golay code, h integral with even coordinate sum}, and its translate by the all-one vector. Closure follows because the binary Golay code is doubly even and self-orthogonal. The even part is integral and has no roots; the odd part has squared norm at least 3 because every physical coordinate is odd. The extended lattice is integral, unimodular and odd with squared minimum 3, so it is the odd Leech lattice; compare the [catalogue coordinate basis](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/oddLeech.html). The determinant and generator containment of the supplied basis are checked exactly.

For a point written as g·1+2c+4h, g∈{0,1}, the class modulo 2O is specified by g, c modulo the all-one Golay word, the binary Golay syndrome of h, and the sum of h modulo 4. These invariants are unchanged by 4c'+8h' in the doubled even sublattice, and by doubling the all-one vector: that complements c and increments h by c, leaving both syndrome and sum unchanged. Conversely, if these four invariants agree, after optional complementation the remaining h difference is a binary Golay word plus twice an integer vector of even coordinate sum, so the difference belongs to 2O.

For each g,c, dynamic programming over the syndrome and sum residue takes all h-coordinate residues modulo 4, minimizing the sum of squares in physical coordinates. The cost at a coordinate with residue a is the square of the least absolute integer congruent to g+2c_i+4a modulo 16. This is exhaustive because changing a coordinate h by 4 leaves all invariants unchanged. Both choices of c up to complementation are taken. The squared depth is the resulting integer cost divided by 16. The program identifies every one of the 2^24 binary coefficient words of the lattice basis bijectively with its invariant class, and applies the exact integer Walsh transform to the two nontrivial functions. No numerical optimization is used in replay.

Run python3 verify.py with SymPy, a C++17 compiler, and about 400 MB of memory.
