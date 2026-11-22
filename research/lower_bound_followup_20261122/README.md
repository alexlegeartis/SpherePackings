# Lower-bound search follow-up in dimensions 19–32

This pass found no strict density improvement over the benchmarks in the [lattice catalog](../../LATTICE_LOWER_BOUNDS_1_48.md) and the [general catalog](../../LOWER_BOUNDS_1_96.md). It gives two exact results that constrain further searches and records the scope of the numerical experiments.

## Wrapped projection of Leech sections

For the 19D, 20D, 21D and 23D antipode normal Grams, the search evaluates every coset of the period in the projected Leech lattice. Choosing arbitrary cosets permits different representatives of each difference in the normal coordinates. The search checks all squared separations and all subsets of the finite quotient through a maximum-clique calculation. It also checks additive subgroups for lattice packings. These four choices of section come from Chen, Hu, Li, Wang and Wu. [Construction and target densities](https://arxiv.org/html/2505.02394v1#S3)

| Dimension | Normal determinant | Highest periodic center density in this family | Highest lattice center density in this family |
|---:|---:|---:|---:|
| 19 | 300 | 0.0889581021 | 0.0577350269 |
| 20 | 125 | 0.1559337389 | 0.0894427191 |
| 21 | 36 | 0.2100363836 | 0.1666666667 |
| 23 | 6 | 0.5004897701 | 0.4082482905 |

For a primitive normal lattice with Gram matrix `K`, the coset labels are the integer pairings `p` modulo `K Z^k`. A representative has squared normal length `p^T K^-1 p`. Each horizontal vector in that coset lifts to a Leech vector with that representative. Its projected squared length is therefore bounded below by `4 - p^T K^-1 p`, or by `6 - p^T K^-1 p` when the norm-4 shell has no vector of those pairings. Conversely each ambient vector gives an upper witness for its coset separation. The two bounds agree for every coset in these four examples. One 19D class requires the explicit norm-6 witness in [leech_normals.json](leech_normals.json), which gives exact projected squared separation 3. The other classes have norm-4 witnesses.

The ambient coordinates are the standard integer Leech coordinates divided by `sqrt(8)`. With `C` the extended binary Golay code, the even part is `4D_24 + 2C`. The odd part is its translate by `(-3,1,...,1)`. The two-coordinate even vectors and the signed octads give the full even norm-4 shell. An odd vector of norm 4 has one coordinate of magnitude 3 and 23 of magnitude 1, with signs prescribed by the Golay code. The code enumerates those vectors exactly and verifies that each normal lattice is primitive by checking surjectivity of the pairing map. [Standard Leech coordinates and short vectors](https://math.berkeley.edu/~reb/papers/leech/leech.pdf)

The exact results cover unions of translates of the displayed periods by points of the projected Leech lattice. The 19D packing has two separate triangles of active site pairs. A continuous search with seven or eight independently moving sites over its fixed period reached squared separations about 2.86667, below the respective requirements 3.41087 and 3.36326. Those runs used floating closest-vector searches and local optimization.

## Quaternary code searches

The first finite batch checked 150 literal shortenings of the original 30D quaternary lattice, 82 signed coordinate identifications, and 346 five-generator subcodes from heuristic candidates. Each used exact distances on the half-grid. The translate and replacement searches yielded no density increase in 27–29D. Re-expressing the three 30D ties in quaternary coordinates also gave no half-grid phase with horizontal squared clearance at least 13/4.

In the binary residue search, 120 bounded annealing runs made about 4.60 billion proposed updates in dimensions 27–32. The closest 28D candidate has 31 weight-12 words, 31 weight-16 words, and the all-ones word. For each weight-12 word, the parity columns on its support span a hyperplane in the six-bit syndrome space. The lifted upper syndrome has to lie outside that hyperplane or the word admits squared separation 3. Seven such equations sum to `0 = 1` for additive lifts. A separate set of 17 equations also sums to `0 = 1` when each of the 64 sites can have an independent lift. The [exact replay](check_lift.py) verifies both contradictions. This rules out the displayed residue set with this period and target minimum, including independent site lifts over all 64 residues.

The numerical metric tests over the three 30D ties and the 26–28D Bacher sections found no first-order density descent. The 26D section has the 117,936-vector shortest shell of T26. Three bounded Voronoi searches there reached squared clearance 6 in the minimum-8 scaling, the value attained by the existing Bacher stack. These numerical searches do not give global covering or local optimality certificates. [T26 catalogue entry](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb4.html)

## Replay

From this directory, run `python3 check_lift.py` and `python3 projected_cosets.py`. The projection replay needs NumPy and SymPy. It compares exact rational squared center densities with [projected_coset_results.json](projected_coset_results.json). Use `--write` to regenerate that file.

## Follow-up to the 29D unrestricted improvement

The newly published [alternating Vardy layers](../vardy28_stack29_20261122/) raise the 29D unrestricted center density to 1. Both verifiers in that package pass. The 29D lattice benchmark stays at `1/sqrt(2)`.

Replacing the vertical line by a triangular lattice of squared minimum 8 would give a 30D center density of `2/sqrt(3)` if the even horizontal phase differences at all three edges of a triangle had Hamming coset depth at least six. For shifts of the form `2b` in either horizontal stock, every same-stock closest pair with nonzero binary weight can realize squared clearance four times that depth. A coordinate with difference 2 modulo 4 permits both signs at the same cost, so the parity condition in `4D_28` does not increase the distance. The three edge syndromes add to zero. In the published code and five other tested circulants, the 78 depth-six syndromes contain no three whose sum is zero, including in each dual code. This also obstructs the direct 31D face-centered cubic extension at the same squared vertical minimum because its shortest edges contain triangles. Run `python3 check_vardy_layers.py` to reproduce the finite checks.

An additional screen retained 1,764 rank-one modifications of the six binary generator matrices with minimum at least eight in both the code and its dual. None provided a joint depth-seven phase or three pairwise depth-six differences. The obstruction applies to the tested horizontal stocks and binary phases. It leaves continuous horizontal phases and other vertical periods open.
