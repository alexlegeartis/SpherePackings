# Deep ternary phases of a norm-six Leech interface in dimension 27

**Result (scope).** Normalize integer Leech coordinates to minimum squared norm 32. Let Λ be the Leech lattice, choose any mark m of usual norm squared 6 (raw squared norm 48), and let H be the kernel of its usual inner product modulo 3. In ((1/3)H)/H there are exactly **six** nonzero phases with squared distance at least 26 from H. They lie on three of the four lines of a single two-dimensional plane over F₃. Consequently **no four phases** have all six pairwise squared distances at least 26. This rules out the fourth vertex needed by the previous [equilateral-interface target](../leech_ternary_interfaces27_20261125/README.md). It does not rule out skew vertical metrics, nontorsion shifts, index-three carriers of all mark types, or general 27D packings. No density record changes.

## Finite reduction

For a lift w/3+a with w∈H, a∈Λ, its projection to the Leech torus is labelled by w mod 3Λ. For each such class make the list S(w) of y∈Λ for which ||w−3y||²<234, in raw integer coordinates. A given lift is deep exactly when the corresponding residue of (m,y) mod 3 is absent from S(w). The sign and an overall additive color offset depend on the representative a; the number of missing colors does not. Replacing w by w+3t translates S by t. The inner product in raw coordinates differs from the usual one by a factor of 8, invertible modulo 3.

We use the ten Co₀ orbits on Λ/3Λ. The minimum *usual* squared norm and number of shortest representatives distinguish them. The source is Conway–Sloane, *Sphere Packings, Lattices and Groups*, Chapter 10, §3.3; the [mmgroup implementation and orbit table](https://mmgroup.readthedocs.io/en/latest/c_interface.html#c-interface-for-file-gen-leech-reduce-mod3-c) explicitly list the ten types. The code establishes the identifying invariants in our coordinates rather than transferring hexadecimal masks between Golay conventions. The search therefore includes all relative positions of w and m: fix each residue orbit in turn and vary the mark over its **whole** shell.

| Minimum, multiplicity | Size of S | F₃ affine rank | Norm-six signed marks admitting at least one deep lift |
|:--|--:|--:|--:|
| (0,1) |1|0|all|
| (4,1) |4602|24|0|
| (6,1) |554|23|2 (parallel to w modulo 3)|
| (8,1) |2095|24|0|
| (10,1)|4603|24|0|
| (12,1)|808|24|0|
| (12,3)|894|24|0|
| (14,2)|1907|24|0|
| (16,9)|4365|24|0|
| (18,36)|765|24|0|

Only test marks orthogonal to w modulo 3, since otherwise w is not in H. For each mark the direct verifier subtracts one point of S from all the others and literally checks which of the three colors occur. For the exceptional orbit there are two missing colors per admissible mark. Rank 24 by itself would rule out *two* missing colors but not one; the explicit shell check in `direct.cpp` is essential to the stronger conclusion.

For reproducibility write each raw Leech point as p·1+2c+4k for p∈{0,1}, the parity-extended cyclic Golay code, and ∑kᵢ≡p mod 2. The norm-48 raw even points have eight entries ±2 and one ±4 outside a Golay octad, or twelve entries ±2 on a Golay dodecad, with the sign parity fixed as in the code. (Three ±4 entries fail the global parity.) The odd points have three entries ±3 or one entry ±5 among the ±1 entries; the residues modulo four fix all signs. These exhaust the possible sums of coordinate squares; there are 16,773,120 points. The S(w) lists use coordinatewise enumeration with a lower bound from the remaining coordinates; the radius in each coordinate is less than 6, so the tested interval about the nearest translate covers every possible lift. `replay.cpp` checks the orbit identifiers a second way using a parity-distance dynamic program and multiplicity counts.

For fixed m the only potentially nonempty projections of deep phases are thus 0 and ±m mod 3Λ. Their inverse image in ((1/3)H)/H is a plane of order nine. The previous note proves that the three squared depths above m are 16/3, 80/3, 32; the opposite fibre is its negative. Above zero the values are 0, 32, 32. The zero and the two shallow opposite points form one subgroup of order three. Thus the graph of allowed differences inside the plane is the complete tripartite graph K₃,₃,₃. Taking one point of any proposed clique as origin proves the stated cap of three everywhere.

## Replay

Compile each C++17 file with optimization and run it: `c++ -O3 -std=c++17 replay.cpp -o replay && ./replay`, then `c++ -O3 -std=c++17 direct.cpp -o direct && ./direct`. The optional target flag `-mpopcnt` speeds up the mark enumeration on x86 supporting it. No external run-time dependencies or input are required.
