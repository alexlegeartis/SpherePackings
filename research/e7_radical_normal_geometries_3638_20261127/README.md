# Full radical E7 layers in dimensions 36–38: exact sharp cap and alternate phase geometries

These are constructions **matching** the maintained densities in dimensions 37 and 38, and an upper bound **only for the stated full-layer carrier**, not an unrestricted upper bound or a new record. The period-free bound also matches the current 36D lower bound. The phase stocks in 37/38 use different normal geometries and their examples here are nonlattices.

## Exact horizontal carrier and phase distances

Use the verified five-block Hermitian carrier [E7^5 syndrome lattice, November 23, 2026](https://github.com/alexlegeartis/SpherePackings/tree/main/research/e7_hermitian35_20261123). Write E for the E7 lattice with squared minimum 2 and determinant of Gram matrix 2. The labels E/2E are (x,e) in F64 × F2 via the included binary-linear parity map, and the norm data modulo 4 are recomputed from the literal Gram matrix. The field uses X^6+X+1. Define the syndrome of a vector in E^5 by

    S = (sum x_i, sum alpha_i x_i, sum e_i), alpha=(0,1,14,15,22).

The tag is packed as s0+64s1+4096e. The map is onto: any two distinct alpha solve the two field equations, and one radical bit supplies the remaining check. Its kernel H has covolume 8192·2^(5/2)=2^(31/2). The 35D carrier H has squared minimum 8 and center density 4.

Let r be the short E7 representative with coefficient mask 82, corresponding to label (0,1). A radical half-phase is any translate of H by

    (b_i r/2)_i + w,  b in F2^5, w integral in E^5.

The squared distance from this translation to H is defined by minimizing over all lattice translates. It depends only on b and the arbitrary integral syndrome S(w). The maximum over syndromes, for every individual b of a given weight, is:

| Weight of b | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|---:|---:|
| Maximum exact squared depth | 6 | 11/2 | 5 | 13/2 | 6 | 15/2 |

For weight 5, **every one of 8192 syndromes** has depth 15/2, not just the maximizer. For any two radical half-phases, subtracting their representatives gives a radical half-phase with mask the XOR of the original masks (and an integral carry). Thus the depth of every even-parity mask difference is at most 6 regardless of the integral carry.

### Replay and completeness

Run `node verify.js`. The exact local table consists of the minima of the E7 Gram norm in all 4^7 coefficient residue classes modulo 4. A positive diagonal factorization of 60 G permits complete integer ellipsoid enumeration through norm 24 (all classes are filled). For fractional mask b and local integral label (x,e), the local cost is the mod-4 minimum for the doubled representative b·r+2·M[x,e]; it is four times the physical squared distance. The exact global DP is an XOR min-plus convolution over the 8192 syndromes. At each site all 128 parity labels are admitted, so this does not omit any integer translate. A paired DP excludes the all-zero local label assignment to verify the horizontal minimum; in the all-zero class a nonzero vector of (2E)^5 has squared norm at least 8, attained by a doubled root. The minimum and all 32 masks are checked independently within the same literal replay.

## Sharp bound for arbitrary normal arrangements

Keep squared sphere separation at least 8. Allow **any** collection of full translates of H above arbitrary points in k-dimensional normal space, k=1,2,3, with every horizontal phase an arbitrary radical half-phase as defined above. In particular the normal arrangement need not be a lattice, and there is no restriction on the period or its shape. This assertion does not cover partial horizontal layers, other fractional alphabets, or a change of radius relative to the fixed horizontal carrier.

Color a layer by parity of b. Two layers of the same color have horizontal squared distance at most 6, so their normal points must be separated by at least √2. Thus each color is a k-dimensional packing with minimum distance √2, whether or not the whole arrangement is periodic. The maximum intensities in dimensions 1–3 are respectively 1/√2, 1/√3 and 1/2 (equal intervals, triangular disks, and Kepler's sphere-packing bound). Adding the two color intensities and converting to center density gives

    delta_(35+k) <= 4 (√2)^k × 2 × I_k
                 = 8, 16/√3, 8√2 in k=1,2,3.

In asymptotic density use the upper bound for each color in large normal boxes, and the uniform full-lattice point-density limit within each horizontal slice; the horizontal lattice shifts have uniform boundary error. This also proves the claim for upper box density. The 1D equality already follows from the [existing order-four radical lift](https://github.com/alexlegeartis/SpherePackings/tree/main/research/e7_hermitian35_20261123).

## New matching 37D and 38D witnesses

The integral syndrome tags 79, 418, 628 have squared depth 6 each and on all three pairwise XORs. Their triple XOR has squared depth 4. Obtain an integral lift of any tag (a,b,e) using labels (a+b,e) and (b,0) in the first two blocks (alpha=0,1), then zeros. In the formulas below the integral shift is its lift; any lift gives the same full layer.

For 37D let P have basis u,v with Gram [[2,1],[1,2]], and put one horizontal layer at every point i u+j v with tag (i mod2)79 XOR (j mod2)418. Zero tag corresponds to even i,j; the next squared normal length there is 8. At normal lengths squared 2 or 6 the nonzero phase has depth 6. Duplicate this entire stock with normal offset (u+v)/3 and horizontal offset (r/2)^5. The normal cross-bank squared distance is at least 2/3, and the horizontal depth exactly 15/2. Thus squared separation is at least 8; the horizontal minimum attains it. The normal area is √3 per point in either bank, so combining the two gives center density 16/√3.

For 38D let P=D3 with basis u=(1,1,0), v=(1,0,1), w=(0,1,1), and tag a site by XOR of 79,418,628 according to the parity of its basis coefficients. Squared normal length 2 gives a single or pairwise tag; the triple-parity class has squared normal minimum 4. The zero tag belongs to 2D3 and cannot give a nonzero normal squared length below 8. The same-bank squared separation is therefore at least 8. Duplicate this stock with normal offset (3/4,1/4,1/4) and the same horizontal all-five half-radical offset. Distances from that normal offset to D3 have squared minimum 11/16: the closest integer point has odd parity and costs 3/16, while changing one coordinate to the other parity costs another 1/2. Adding the uniform horizontal 15/2 gives more than 8. The normal covolume of either stock is 2; combining both gives center density 8√2.

In both dimensions the doubled normal offset is neither in P nor in the translated stock, proving nonclosure under doubling for a packing containing zero. Distances vary continuously with the normal offset, and the cross-stock inequalities are strict, so the equality-density examples include an open neighborhood of offsets without a spacing loss. This slack is a useful starting point for a genuinely different phase type; it cannot give a density gain while all layers remain radical half-phases.

Benchmarks: [current SpherePackings table](https://github.com/alexlegeartis/SpherePackings/blob/main/LOWER_BOUNDS_1_96.md), [37D primitive section](https://github.com/alexlegeartis/SpherePackings/tree/main/research/e7_37_section_20261123), [38D coupled E7 lattice](https://github.com/alexlegeartis/SpherePackings/tree/main/research/e7_lattice38_20261123), and [Cohn’s table](https://cohn.mit.edu/sphere-packing/#data). The tied constructions do not improve those entries.

