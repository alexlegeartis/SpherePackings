# An 11D 72-site packing with a nonbinary phase

**Density tie, not an improved bound.** Here is a 72-site motif with period 2Z^11, minimum squared separation 4, and center density 9/256. It is not isometric to any packing on a translated integer grid: one squared interpoint distance is an integer plus 1/4. We make no priority claim for free-site perturbations of this classical family.

## Construction

Start with the normalized order-12 Paley Hadamard matrix. Its first row is all +1. For the other rows indexed by i in F11, the first column is +1, the diagonal of the remaining block is −1, and an off-diagonal entry is +1 exactly when i−j is a nonzero square. For each pair of distinct rows take their coordinatewise product and its negative. Retain the 66 products with first coordinate +1 and delete that coordinate. Take each −1 to 1 and +1 to 0. The result is a 66-word binary core of weight six.

Add five binary words of weight two on the disjoint pairs encoded as 1028, 520, 272, 160, 66 (least significant bit is coordinate zero). Add one more word based on the all-ones vector, shifting coordinate 10 by +1/2, modulo the period 2 (equivalently take that coordinate to −1/2). The exact minimum of the wrapped squared distances between distinct motifs is 4, and within one motif the period lattice has minimum 4. Thus all periodic images are accounted for and the center density for balls of radius one is 72/2048.

Run `python3 verify_tie.py`. The finite calculation verifies the matrix, all pair distances including the noninteger witness, and the density. For the historical binary benchmark see [Conway and Sloane, Quaternary constructions for the binary single-error-correcting codes of Julin, Best and others](https://neilsloane.com/doc/Me185.pdf), and [Nebe–Sloane's 11D table](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/density.html).

## What can and cannot be changed locally

The nonbinary phase is a genuine change of the motif, not evidence that the 66-site core alone is deformable. In fact the core is strictly locally maximal at its own density at fixed motif count and chosen period, even with joint metric and site moves. This does not constrain a larger supercell, core replacement, or a separated phase such as the one above.

To see this, use a variable Gram metric at the cubic cell and normalize the minimum squared distance to 4. Write S for its derivative and u_i for the displacement of site i. The native period contacts 2e_j force S_jj≥0. Nonincreasing volume forces tr S≤0, so all diagonal entries vanish. A pair i,j of core words differing in exactly four coordinates has all 16 nearest displacements v, with independent choices of signs on its four supports. The derivatives v^T S v +2 v·(u_i−u_j) are nonnegative and have average zero. They are therefore all zero. Mixed sign differences force S_ab=0 whenever {a,b} occurs in any support, which the finite check verifies for every pair of coordinates. Individual sign differences then force the coordinate-a motion to agree along every core contact whose support contains a. Those 11 contact graphs are connected, so all motions are common translations. Fix one site to quotient translations. A contradicting sequence of distinct nearby admissible packings of no larger volume would normalize, by compactness, to a nonzero direction of the ruled-out kind.
