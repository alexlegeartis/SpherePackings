# Two ternary Leech interfaces in dimension 27, and a different interface with a deep triangle

**Scope.** In the coordinate realization of the Leech lattice below, take either index-three sublattice H given by one marked coordinate y23 = 0 modulo 3 or by y22 + y23 = 0 modulo 3. For either carrier, any rank-three vertical lattice V, and any homomorphism V → (H/3)/H (here H/3 means (1/3)H), stacking a copy of H at each site of V according to that phase map gives 27D center density at most **9/16**. This includes arbitrary vertical metrics and non-surjective phase maps. No comparable statement is asserted for all index-three carriers, for nonlinear assignments, or for nonternary phases. Indeed the last section exhibits three compatible deep phases for a different carrier.

The family ceiling 9/16 < 1/√3 is below [Bacher's 27D lattice](https://doi.org/10.1007/s002220050180). [Vardy's unrestricted 27D packing](https://doi.org/10.1090/S0002-9947-99-02169-8) attains 1/√2. **No record changes.**

## Coordinates and complete finite test

Use the usual integer coordinates for the Leech lattice Λ, with minimum squared norm 32, covolume 8^12 and center density one. With G the extended binary Golay code, all points are
`y = p 1 + 2c + 4k`, p = 0 or 1, c ∈ G in 0/1 representatives, and the sum of k is p modulo 2. The twelve generator rows of G are the shifts of the cyclic binary polynomial with mask 2787, parity-extended at position 23. The constraints defining H are nontrivial since 8e23 ∈ Λ. Both H contain a vector 4(e0+e1) of squared norm 32.

Put u=8e23 modulo H, generating Λ/H. The 3^24 horizontal phases form fibres of three under projection to the Λ torus; the three phases in each fibre differ by u. Write D for distance squared from H. The following is the complete count of fibres having ANY member with 9D > 192; all remaining fibres have 9D ≤ 192 throughout:

| Sorted triple of 9D | one-coordinate mark | two-coordinate mark |
|---|---:|---:|
| (0,288,288) | 1 | 1 |
| (32,224,224) | 69,876 | 93,150 |

The unique first fibre is {0,u,−u}. Explicitly write a phase as 8a/3 + tu, t ∈ {−1,0,1}. On the unmarked coordinates a is an arbitrary vector of trits represented by 0, ±1. On the single mark a23=0; on two marks take (a22,a23)=(r,−r) with r=0,±1. These are distinct: if two differed by H, the numerator difference would be integral and hence all a coordinates equal mod 3, then their t coordinates also agree. They lie in (1/3)H and give exactly 3^24 labels.

The verifier partitions support masks under explicitly tested Golay code coordinate permutations fixing the marks. For each support, binary sign flips from Golay codewords zero on the marks are valid Leech isometries: the sign flip changes the total-coordinate sum by a multiple of 8 (the Golay code is doubly even and self-dual), and preserves the Golay residue condition. On that support their projections form a binary linear subspace; enumerate complementary sign masks by Gaussian elimination and multiply by its size. No reliance on transitivity claims about Mathieu groups: the program generates every support orbit literally. For each label and word c,p, the unmarked coordinates must lie in p+2c+4Z. A marked coordinate also lies in 3Z; for two marks impose zero total mod 3 jointly. For each factor retain the nearest cost in each carry parity, then impose the single global parity by changing the cheapest parity if necessary. Calculations use squared integer residuals with denominator 9 cleared. A second distance procedure directly propagates carry parity (and the sum modulo three for two marks) on a sample of labels. The small integer search ranges in the direct check cover nearest values in each residue class: unmarked k has period 2 and marked k has period 6; the marked coordinate targets have absolute value at most 32/3. The combined checks count 282,429,536,481 labels in each case.

## From fibres to arbitrary metrics

Let b be the squared length of a primitive shortest vector v in V. Since 3v has zero phase, b ≥ 32/9. If b ≥ 64/9, the classical sharp 3D Hermite inequality gives covol(V)^2 ≥ b^3/2 ≥ 131072/729.

Otherwise the phase of v must be ±u: every other phase has D ≤ 224/9 and needs vertical squared length at least 64/9. Orthogonally project V along v to a rank-two lattice. For each nonzero projected vector of squared length t, the lifts above it have horizontal phases cycling through one of the triples above. If all three depths are at most 192/9, use a lift whose longitudinal displacement has magnitude at most |v|/2. Then t ≥ 32 − 192/9 − b/4 > 80/9. Otherwise one depth in the triple is at most 32/9. In that particular phase, choose a lift with longitudinal displacement at most 3|v|/2. Then t ≥ 32 − 32/9 − 9b/4 > 112/9. Thus every nonzero projected vector has squared length at least 80/9. The planar Hermite inequality gives
`covol(V)^2 ≥ b (3/4)(80/9)^2 ≥ 153600/729`, even larger. Both cases imply covol(V) ≥ 256√2/27. The unchanged horizontal radius √8 and horizontal index 3 then give center density ≤ (16√2/3)/(256√2/27)=9/16.

## Constructive contrast: a norm-six mark has an independent triangle

For this paragraph only use the usual normalization of Λ with squared minimum four (divide the preceding integer coordinates by √8). Choose m ∈ Λ with squared norm six and set H={y∈Λ : (m,y)=0 mod 3}. In the fibre m/3+Λ modulo H the three squared depths, **after conversion back to the preceding scale**, are 16/3, 80/3 and 32. The last two meet the required depth 26 for nearest interfaces of the proposed √3 D3 stack. Together with u from Λ/H this gives a triangle of phases of affine dimension two over F3: choose a lift z on one of the two deep lifts, then choose the sign of u so z−u lies on the other. All three edge depths of {0,u,z} are at least 80/3. It does not supply the independent fourth vertex or a corrected record.

Here is a short exact verification. For y∈Λ the numerator n=||m−3y||² is 6−6(m,y)+9||y||²; in each H-coset n mod 18 is fixed. Neither n=12 nor n=18 is possible: for ||y||² ≥ 6 the reverse triangle inequality gives n ≥ 24, and for ||y||²=4 either value would give ||m−y||²=0 or 2, impossible. A norm-32 integer-coordinate raw root with normalized scalar product (m,y)=2 gives n=30, and one with scalar product 1 gives n=36. For explicit witnesses take m with a 4 at coordinate 0, entries ±2 on a Golay octad avoiding 0 (one negative), zeros elsewhere, in the raw coordinates. Roots 4(e0+ej), j outside the mark, and 4(e0−ei), i with a +2, give scalar products 2 and 1. The value n=6 occurs at zero. The three residues 6,12,0 mod 18 are distinct and attain the stated minima. The nonzero Λ/H cosets also have squared depth 32 in the original scale (one of the displayed roots is nonzero mod H, and its negative supplies the other). The three phases have order dividing three since (m,m)=6; z does not belong to Λ/H since m/3 is not a lattice vector.

## Replay

`c++ -O3 -std=c++17 census.cpp -o census && ./census` on Linux with GCC/Clang and the standard C++ headers as supplied by libstdc++. No external datasets or floating-point computations are used. No generated tables need to be downloaded. The source verifies both fibre counts and sample distances before completing.
