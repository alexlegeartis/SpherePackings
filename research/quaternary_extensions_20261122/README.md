# Lattice and periodic packing searches in 29–31 dimensions

This first pass uses the repository’s [1–48D lattice catalog](../../LATTICE_LOWER_BOUNDS_1_48.md) and [general lower-bound catalog](../../LOWER_BOUNDS_1_96.md) as comparison inputs. The existing 29D and 30D quaternary lattices already attain the general benchmarks. In 31D both benchmarks have center density `3^15 / 2^(47/2)`, about 1.2095224. We used these three dimensions to test whether adding centers to a fixed period or stacking a new layer could improve both benchmarks.

## Results

Three 30D lattices have center density **1** and kissing numbers **123,212**, **124,748**, and **125,260**, respectively. The original quaternary lattice has the same density and kissing number **128,844**. The four different contact counts prove that the four lattices are pairwise nonisometric. The corresponding occupied-volume density is at least **2.19153534478302e-5**. The best-density entries in both catalogs still stand after this pass.

| Search | Cases | Result |
|---|---:|---|
| Add one translate at squared separation 4 | 139 | No admissible translate on the half-grids or tested one-axis quarter-grids |
| Replace one coset of an index-two sublattice at squared separation 4 | 252 | Admissible individual translates, zero mutually compatible pairs on the half-grids |
| Add one translate at squared separation 31/8 | 1,772 | No admissible translate on the tested one- and two-axis quarter-grids |
| Stack a lattice over the 29D or 30D seed with a half-grid phase | 4 | Twelve phases tie the 30D density, zero 31D extensions at either height |
| Stack using one- or two-axis quarter-grid phases | 900 | No admissible phase at squared height 7/16 over 29D or 5/8 over 30D |

The fixed-period tests start with the [quaternary witnesses](../quaternary_29_30_20261122/witness.json). The 29D test folds equal coordinates 2 and 3 into one coordinate of squared-length weight two. The 31D tests keep the six original binary residue generators and restore either omitted coordinate. Of the 64 choices of high bits on a restored coordinate, eight preserve the minimum for each choice of coordinate. The unshifted translate tests include all sixteen resulting seeds. The shifted and index-two tests use the two 31D seeds with zero high bits on the restored coordinate. The index-two test covers all 63 binary residue hyperplanes for each of its four seeds. Its graph asks whether the retained coset and two added cosets can coexist. Three cosets would multiply the corresponding full seed density by `3/2`, enough to beat the listed benchmark in each dimension. The smaller-radius test also has a strict benchmark gain if feasible. Two cosets multiply the seed density by `2(31/32)^(n/2)`.

## Lattice witnesses and verification

Write a horizontal phase as `(x + 2z)/2`, with `x` and `z` binary masks. The folded coordinate uses the same scale in each of the two original axes. The six-bit syndrome of `z` suffices because changing `z` by a torsion word leaves the coset unchanged. To reconstruct a phase from the report, take any binary `z` with the recorded syndrome. The following three phases use a squared layer height of `1/2` over the 29D lattice.

| Low-bit mask `x` | High-bit syndrome | Kissing number |
|---:|---:|---:|
| 66,121,728 | 12 | 123,212 |
| 217,116,672 | 57 | 124,748 |
| 66,121,728 | 8 | 125,260 |

The horizontal squared distances of the phase and its double from the 29D lattice are exactly `7/2` and `2`. Adding the vertical squared lengths `1/2` and `2` gives 4 in each case. Every layer of absolute index at least three has vertical squared length at least `9/2`. The horizontal layer itself has minimum 4 and covolume `sqrt(2)`, so the 30D lattice has minimum 4, covolume 1 and center density 1. The verifier counts all contacts in layers 0, ±1 and ±2 by summing over the literal integer representatives of each coordinate. It verifies the original 30D contact count both in its quaternary representation and in its layered representation.

For an unshifted half-grid phase, the first horizontal squared distance is a multiple of `1/4`. The second is an integer because doubling makes the binary residue zero and all seed residue weights are divisible by four. Four times the phase lies in the horizontal lattice, so the squared height is at least `1/4`. The third layer has the first layer's horizontal distance. Consequently the least feasible squared height is a multiple of `1/4`. Testing `1/4` and `1/2` suffices for beating the benchmarks using these half-grid stacks. The quarter-grid tests check only the displayed heights and offset patterns.

## Exact search and replay

For a fixed offset pattern, clearing the binary residue pivot bits leaves one low mask per quotient class. Each of its 64 high syndromes gives one phase modulo the seed lattice. The C++ search computes the attainable syndromes below the distance threshold by a weighted subset-sum dynamic program and removes them. The remaining choices are intersected across all seed residues. For the replacement search, a difference of two pivot-cleared low masks remains pivot-cleared. The phase pair is therefore compatible exactly when the difference, including the subtraction borrow in the high syndrome, is in the same admissible set. The code tests whether twice a phase lies in the period when checking index-two lattice extensions. For the two-coset replacement test it also checks the cyclic order-four and the two-generator index-four conditions. The Python replay recomputes witness distances by a separate integer shortest-path recurrence.

From this directory, use Python 3.10+ and a C++20 compiler. Each command compares its output with the corresponding checked-in JSON. Add `--write` to regenerate it. `--jobs` controls parallelism for the searches that accept it.

```sh
python3 replay.py --jobs 16
python3 replay.py --repack --jobs 16
python3 replay.py --near-contact --jobs 16
python3 stack.py
python3 stack.py --quarter-grid --jobs 16
```

[Translate results](results.json) · [Replacement results](replacement_results.json) · [Smaller-radius results](near_contact_results.json) · [Half-grid stacks](stack_results.json) · [Quarter-grid stacks](quarter_grid_stack_results.json)

## Candidates for the wider 1–48D search

The lattice-only and general baselines differ most in 28D (50.0%), 45D (41.0%), 46D (25.9%), and 20D (24.7%), measured as the gain needed to reach the general density from the current lattice density. The two catalogs also give separate targets for 27D, 43D, and 44D. On the lattice side, 27–28D support trying different residue codes or torsion codes, and 44–46D support testing additive gluing over the known 48D sections. For general packings, the 30D ties give distinct contact geometries to test under continuous deformations. Farther up the range, the record-setting 38–42D sections and the six-translate 43D packing are explicit seeds. [Dorofeev, Sun and Wang, Theorem 1](https://arxiv.org/html/2607.20359v4#Thmtheorem1)

The general problem is solved in 1–3D, 8D, and 24D. The lattice problem is also solved through 9D, so dimensions 4–7 and 9 are targets for the general search. [Dutour Sikirić and van Woerden](https://arxiv.org/abs/2508.20719) · [Exact optima and sources](https://arxiv.org/html/2607.20359v4#S1)
