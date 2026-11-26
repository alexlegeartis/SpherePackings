# 13D: exact limits of integral planar OBK fibers and real glides of their unit rows

SpherePackings round 94, November 26, 2026. Public base for the comparison: center density 9/256, occupied density 0.0320142921603497... in [Cohn’s table](https://cohn.mit.edu/sphere-packing/#data) and [the project’s 13D entry](https://github.com/alexlegeartis/SpherePackings/blob/22ed04d3ae7c517367707742dc9fd6d45cf9daeb/LOWER_BOUNDS_1_96.md). Minimum distance 2 throughout. The objective was a 13D improvement from rhombic planar layering. We obtain exact obstructions and an equality realization, not a new density record. Neither statement classifies general 13D packings.

## 1. Integral phases: a global result without periodicity or a catalogue hypothesis

Let Y be any set of layer sites in the Euclidean plane. Above each y put a fiber C_y+2Z^11, with C_y ANY 72-element binary length-11 code of minimum Hamming distance 4. Choices at different sites may be unrelated. Integral shifts and signed coordinate permutations are absorbed into C_y; all such fibers use the same coordinate axes and lattice 2Z^11. Assume the resulting 13D separation is at least 2. The upper planar point density is at most 1, so the center density of the 13D packing is at most 72/2048=9/256.

First, |y−z|≥1: otherwise every binary pair across the two codes has Hamming distance ≥4 (including disjointness of the codes). Their union has size 144 and distance 4. Puncturing one bit gives a binary (10,144,≥3) code, contradicting the radius-one Hamming bound 144·11>1024.

Draw an edge when |y−z|<√2. There is no triangle: again its three codes must be disjoint, now with union size 216 and distance ≥3, contrary to 216·12>2048. Edges of this straight-line graph cannot cross: a convex quadrilateral has an angle at least a right angle; the diagonal across that angle therefore has squared length at least 2 since both adjacent sides have length at least 1. A vertex cannot lie in the interior of a nonincident edge of length <√2. Every finite subgraph is planar and triangle-free and has at most 2N edges.

For finitely many sites truncate each site's Voronoi region to its disk of radius R=1/√2. It has area at least π/2 minus c for each neighbor at distance <√2, where c=π/8−1/4 is the largest circular cap cut by a bisector at distance at least 1/2 from the site. Distinct truncated regions have disjoint interiors. Their total area is ≥Nπ/2−2|E|c ≥N. Apply this to the sites in a large square: the truncated regions are in the square enlarged by R in each direction. Taking the limit proves the density bound. Uniform site count ≤(L+2R)^2 in any L-square and uniform horizontal 72/2048 intensity also give the claimed Euclidean upper density without periodicity.

## 2. Real glides: a different theorem, with a fixed 12D unit-row bank

Let C be one of the 96 listed binary (11,72,4) codes of [Östergård–Baicheva–Kolev](https://users.aalto.fi/~pat/72.html). Let v be a binary translation with every pair from C and C xor v at Hamming distance ≥3, so the 11D fibers can alternate at unit spacing along the first exterior direction. A full 12D row is the periodic packing
`P=(B+2Z^12), B={(c,0):c∈C} ∪ {(c xor v,1):c∈C}.`
This is feasible exactly for 27 of the fixed catalogue representatives, always v=all ones and unique. This does NOT claim that all possible unit rows have this form.

Parallel copies of ONE such 12D P are now placed in 13D, allowing a separate arbitrary REAL translation in its entire 12D hyperplane on every row (including arbitrary axial stagger), and arbitrary nonuniform normal gaps. The orientation and the bank P are common, and rows are full. Every gap is at least 1. Therefore center density ≤144/4096=9/256. The proof is substantially stronger than a half-grid search.

Specifically, E=(B xor B) is the binary residue set of P−P, and E+2Z^12 has Euclidean covering radius exactly √3. For each of the 4096 binary corners u, either u∈E or there exist four elements of E for which the supports of u xor e are pairwise disjoint. The standalone finite proof below finds masks of weight ≤3 throughout all 27·4096 corners. To see why this proves the continuum claim, choose a nearest binary corner u coordinatewise to arbitrary real t modulo 2. Put a_i=1−2 dist(t_i,u_i+2Z)∈[0,1]. The squared distance to e+2Z^12 is exactly
`Σ_i (1−a_i)^2/4 + Σ_{i:(u xor e)_i=1} a_i.`
For four disjoint supports one weighted sum is at most Σa_i/4; if u∈E choose the empty support. The resulting upper bound is `3−Σ_i a_i(1−a_i)/4 ≤ 3`. The all-half points attain 3 because every binary word then costs 12/4. Consequently the minimum distance between two shifted copies of P in parallel hyperplanes is at most sqrt(3+h²), implying h≥1. The same argument gives an exact, quantitative loss from maximum for an intermediate phase, not merely local stationarity.

At any maximizer each a_i equals 0 or 1; all coordinates are integer or half-integer. For a fractional support F and complementary parity bits b the squared distance is |F|/4+min_{e∈E} wt((b xor e) outside F). Only |F|=0,4,8,12 can give 3. For EVERY bank, direct projected Hamming-ball checks give respectively 0,0,0,1 possibilities after identifying signs on F. Equivalently all maximizers in the torus are all-half vectors; their independent coordinate signs represent different points modulo 2, but have the same costs. The homogeneous construction P+j(1/2,...,1/2) at normal heights j is valid and attains 9/256. This is an equality witness for (2), not asserted to be inside the integral planar family of (1). No restriction on arbitrary horizontal rotations/changes of bank, partial rows, or arbitrary 13D packings follows from (2).

Independent checks: one native constrained C++ Hamming-distance-transform search found no half-grid triangular winner; the standalone JS calculation below uses direct pair XOR and disjoint-error covers and was separately replayed on all codes. The corner counts, by binary distance to E (0,1,2), are 39924, 55296, 15372. The 27 viable labels are printed by the script.


[Standalone replay](replay.js) verifies the real-glide assertions on the authors’ binary code catalogue; run `node replay.js /path/to/72_4.txt`.
