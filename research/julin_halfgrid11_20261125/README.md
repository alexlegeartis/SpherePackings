# A half-grid extension test for the 11D Steiner packing

Normalize the period to twice the integer lattice and the required distance to 2. For the specific 66 binary Steiner sites below, retain **any three of their six cyclic orbits** (33 specified sites). Even after completely replacing every other site, a packing with this period, these retained sites, and a motif on the half-integer grid has at most **72** sites per period. Equality is attained by the classical construction. This eliminates one direct attempt to beat P11a by increasing the motif to 73. It is not a bound for arbitrary real motifs, larger or deformed periods, arbitrary choices of 33 sites, or unrestricted 11-dimensional sphere packing. The Steiner sites are from [Leech and Sloane (1971), §2.6](https://doi.org/10.4153/CJM-1971-081-3), not a new construction.

## Exact reduction

Number coordinates modulo 11 and use the support {0,1,3,4,5,9}. Form the 12 rows consisting of zero and the cyclic shifts of its indicator. The 66 distinct pairwise XORs form six cyclic orbits of size 11. The five disjoint pair indicators on coordinates (0,1), (2,3), …, (8,9), together with the weight-ten indicator of the first ten coordinates, extend the 66 to a binary motif of 72 at squared distance at least 4. Its center density is 72/2¹¹ = 9/256.

Represent a half-grid point by two 11-bit masks h,b, with coordinates bᵢ+hᵢ/2 modulo 2. For points (h,b) and (h′,b′) set H = h XOR h′. The exact squared flat-torus distance multiplied by 4 is the population count of H plus four times the number of positions outside H where b and b′ disagree.

Retaining three orbits leaves 33 fixed sites. For each of the 20 choices the verifier explicitly generates every eligible half-grid point, then asks for a clique of **40** additional points, with adjacency if and only if the distance numerator is at least 16. The candidate pools contain 981, 1,168, or 1,344 points depending on the orbit choice. In every case the clique search proves that 40 additions are impossible. The old 72-point construction supplies 39 additions.

For generation, the signs on half-integral coordinates do not affect any distance to a fixed binary site. Enumerate the half-coordinate support and the binary coordinates on its complement, then restore all allowed signs. The clique verifier greedily partitions each induced candidate set into independent color classes to bound the possible size of an extension, and branches exhaustively when this bound is not sufficient. No search cutoff, floating-point tolerances, or external solver is involved. Run `c++ -std=c++17 -O3 verify.cpp -o verify && ./verify` (without `-DNDEBUG`).

**Revisit condition.** A denser construction must leave at most two entire orbits unchanged, leave the half grid, change the period, or use a different scaffold. A continuous local search did not certify an improvement and is not part of this theorem.
