# Lattice-only sphere-packing lower bounds, 49D–96D

This catalog separates lattice packings from the [unrestricted 1–96D table](LOWER_BOUNDS_1_96.md).
Each row compares a prior published lattice or an explicitly named derived lattice with
the best verified repository lattice found in this review. A section lattice may be
included even when adding extra translates yields a denser **nonlattice** packing;
only the lattice period counts here. Lattice existence proofs are marked as such.

The center density is δ = (m/4)^(n/2) / sqrt(det G), where m is a proved lower
bound on squared minimum and G is a lattice Gram. Ordinary density is Vₙ δ.
Numerical columns are downward-rounded lower bounds. Gains are rounded down and use
the upper end of a cited rounding interval when supplied; a dash means no
verified increase is claimed. Prior values are the strongest comparisons found
in the cited source audit, not a claim of exhaustive historical priority.

| Dim. | Prior lattice center density ≥ | Repository lattice center density ≥ | Gain over prior | Current occupied density ≥ |
|---:|---|---|---:|---:|
| 49 | 8.41705609801e+3 · [Leech--Sloane P48q x classical 1D lattice](https://neilsloane.com/doc/Me103.pdf) (derived) | 1.45787688119e+4 · [P48q half-coset lift](research/covering_and_deep_cosets_20261115/DEEPER_PROOF.md) | +73.20% | 7.15135575565e-9 |
| 50 | 4.85958960397e+3 · [Leech--Sloane P48q x classical 2D lattice](https://neilsloane.com/doc/Me103.pdf) (derived) | 1.45787688119e+4 · [P48q deep plane with fitted vertical metric](research/adaptive_layers_20261115/results.json) | +200.00% | 2.52244663789e-9 |
| 51 | 2.97587872226e+3 · [Leech--Sloane P48q x classical 3D lattice](https://neilsloane.com/doc/Me103.pdf) (derived) | 1.54631194323e+4 · [P48n quaternionic extension](research/p48n_hurwitz_20261122/README.md) | +419.61% | 9.34489307985e-10 |
| 52 | 5.23755238637e+3 · [Leterrier54, Corollary 1.4 (2022) / Blichfeldt section to 52D](https://doi.org/10.1007/s40993-022-00321-0) (derived) | 1.89383762205e+4 · [P48n quaternionic extension](research/p48n_hurwitz_20261122/README.md) | +261.58% | 3.95931817049e-10 |
| 53 | 1.75328256422e+4 · [Leterrier54, Corollary 1.4 (2022) / Blichfeldt section to 53D](https://doi.org/10.1007/s40993-022-00321-0) (derived) | 2.31946791485e+4 · [P48n quaternionic extension](research/p48n_hurwitz_20261122/README.md) | +32.29% | 1.66176494990e-10 |
| 54 | 6.03065591597e+4 · [Leterrier54 (2022), Corollary 1.4](https://doi.org/10.1007/s40993-022-00321-0) | — | — | 1.46699394710e-10 |
| 55 | 9.04598387395e+4 · [Repaired order-three Leterrier lamination](research/tilted55_20261116/PROOF.md#external-source-correction) (derived) | 9.13338671435e+4 · [Cross-coset Leterrier54 tilt](research/cross_coset55_20261116/PROOF.md) | +0.96% (repaired comparison only) | 7.47532795256e-11 |
| 56 | 1.96213213811e+5 · [Shioda MW56](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | — | — | 5.35530073542e-11 |
| 57 | 9.81066069059e+4 · [Shioda MW56 x classical 1D lattice](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) (derived) | 1.79117338809e+5 · [MW56 depth-7 half-coset lift](research/packing_binary57_58_20261117/README.md) | +82.57% | 1.61599884109e-11 |
| 58 | 1.51348908966e+5 · [Classical Gamma72/Leech-block section](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) (derived) | 1.88806253021e+5 · [MW56 depth-7 complex-multiplication lift](research/packing_binary57_58_20261117/README.md) | +24.74% | 5.58243739609e-12 |
| 59 | 2.69064727051e+5 · [Γ72 section normal to √2 κ13](research/covering_and_deep_cosets_20261115/PRIORITY.md) (derived) | 2.18311667762e+5 · [Explicit bare section of Ne64](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#5963d-bare-ne64-sections) | — | 2.58516365378e-12 |
| 60 | 5.39027741949e+5 · [Elkies E60, reported by Vardy (1999), §4](https://doi.org/10.1090/S0002-9947-99-02169-8) (reported floor) | 6.21378370370e+5 · [Γ72 section normal to √2 K12](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#60d-coxetertodd-section) (derived) | +14.48% | 1.92394316159e-12 |
| 61 | 1.52460068632e+6 · [Elkies E61, reported by Vardy (1999), §4](https://doi.org/10.1090/S0002-9947-99-02169-8) (reported floor) | 1.46447918819e+6 · [Explicit bare section of Ne64](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#5963d-bare-ne64-sections) | — | 1.50881995002e-12 |
| 62 | 4.58980044837e+6 · [Elkies E62, reported by Vardy (1999), §4](https://doi.org/10.1090/S0002-9947-99-02169-8) (reported floor) | 4.14217265951e+6 · [Explicit bare section of Ne64](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#5963d-bare-ne64-sections) | — | 1.44018435505e-12 |
| 63 | 1.58722134875e+7 · [Elkies E63, reported by Vardy (1999), §4](https://doi.org/10.1090/S0002-9947-99-02169-8) (reported floor) | 1.24265179785e+7 · [Explicit bare section of Ne64](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#5963d-bare-ne64-sections) | — | 1.56659999939e-12 |
| 64 | 4.30467210000e+7 · [Nebe Ne64](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ64b.html) | — | — | 1.32606329668e-12 |
| 65 | 3.35544320000e+7 · [Classical Gamma72/Leech-block section](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) (derived) | 4.08377052485e+7 · [Ne64 deep order-three lift](research/novel_routes_20261113/higher_parent/ne65_verified.json) | +21.70% | 3.89626482463e-13 |
| 66 | 7.74906413908e+7 · [Γ72 Leech-block section](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) (derived) | — | — | 2.27252865205e-13 |
| 67 | 1.89812531248e+8 · [Γ72 Leech-block section](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) (derived) | — | — | 1.69831045776e-13 |
| 68 | 5.36870912000e+8 · [Γ72 Leech-block section](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) (derived) | — | — | 1.45479212622e-13 |
| 69 | 1.51850024998e+9 · [Γ72 Leech-block section](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) (derived) | — | — | 1.23719412350e-13 |
| 70 | 4.95940104901e+9 · [Γ72 Leech-block section](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) (derived) | — | — | 1.20626332966e-13 |
| 71 | 1.71798691840e+10 · [Γ72 Leech-block section](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) (derived) | — | — | 1.23869490986e-13 |
| 72 | 6.87194767360e+10 · [Nebe Γ72](https://arxiv.org/abs/1008.2862) | — | — | 1.45861310696e-13 |
| 73 | 3.43597383680e+10 · [Nebe Gamma72 x classical 1D lattice](https://arxiv.org/abs/1008.2862) (derived) | 6.87194767360e+10 · [Γ72 deep half-coset lift](research/gamma72_depth6_20261122/README.md) | +100.00% | 4.26462994945e-14 |
| 74 | 1.98376041960e+10 · [Nebe Gamma72 x classical 2D lattice](https://arxiv.org/abs/1008.2862) (derived) | 7.93504167841e+10 · [Γ72 deep half-coset lift](research/gamma72_depth6_20261122/README.md) | +300.00% | 1.43007109517e-14 |
| 75 | 1.21480019999e+10 · [Nebe Gamma72 x classical 3D lattice](https://arxiv.org/abs/1008.2862) (derived) | 9.71840159992e+10 · [Γ72 deep half-coset lift](research/gamma72_depth6_20261122/README.md) | +700.00% | 5.05260043732e-15 |
| 76 | 8.58993459200e+9 · [Nebe Gamma72 x classical 4D lattice](https://arxiv.org/abs/1008.2862) (derived) | 1.37438953472e+11 · [Γ72 deep half-coset lift](research/gamma72_depth6_20261122/README.md) | +1500.00% | 2.04778582360e-15 |
| 77 | 2.42960039998e+10 · [Bachoc–Nebe L80 section normal to 2A3](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#77d-bachocnebe-section) (derived) | 1.58700833568e+11 · [Half-root Γ72–D4 tilt](research/gamma72_root_tilt_20261122/README.md) | +553.19% | 6.73268761846e-16 |
| 78 | 7.93504167841e+10 · [Bachoc--Nebe L80 bare section, normal rank 2](research/cross_dimension_lower_20261115/l80_results.json) (derived) | 2.11601111424e+11 · [Half-root Γ72–D4 tilt](research/gamma72_root_tilt_20261122/README.md) | +166.66% | 2.53967358970e-16 |
| 79 | 2.74877906944e+11 · [Bachoc--Nebe L80 perpendicular to a minimum vector](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf) (derived) | 3.27807610890e+11 · [MW80 primitive section normal to 162/41](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#79d-bare-mw80-section) | +19.25% | 1.10606606099e-16 |
| 80 | 1.21904190160e+12 · [Shioda MW80](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | — | — | 1.14912817562e-16 |
| 81 | 6.09520950801e+11 · [Shioda MW80 x classical 1D lattice](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) (derived) | 1.09951162777e+12 · [Bachoc–Nebe deep half-coset lift](research/packing81_verified_20261117/README.md) | +80.38% | 2.87777416803e-17 |
| 82 | 3.51907085022e+11 · [Shioda MW80 x classical 2D lattice](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) (derived) | 1.26960666854e+12 · [A2 section of the 83D extension](research/l80_a2_section_20261122/README.md) | +260.77% | 9.17033181889e-18 |
| 83 | 2.74877906944e+11 · [BCS88 section normal to √2 D5](research/bcs83_comparator_20261117/README.md) (derived) | 1.55494425598e+12 · [Bachoc–Nebe parity extension](research/packing83_85_deep_verified_20261117/README_DEEP.md) | +465.68% | 3.08087078338e-18 |
| 84 | 5.49755813888e+11 · [BCS88 root-normal section (4 normals)](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) (derived) | 2.19902325555e+12 · [Bachoc–Nebe parity extension](research/packing83_85_deep_verified_20261117/README_DEEP.md) | +300.00% | 1.18808157306e-18 |
| 85 | 1.09951162777e+12 · [BCS88 root-normal section (3 normals)](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) (derived) | 3.10988851197e+12 · [Bachoc–Nebe parity extension](research/packing83_85_deep_verified_20261117/README_DEEP.md) | +182.84% | 4.55474871517e-19 |
| 86 | 2.53921333709e+12 · [BCS88 root-normal section (2 normals)](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) (derived) | 6.61251878039e+12 · [Primitive MW92 kernel section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | +160.41% | 2.61014583751e-19 |
| 87 | 6.21977702395e+12 · [BCS88 root-normal section (1 normals)](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) (derived) | 1.93579136475e+13 · [Primitive MW92 kernel section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | +211.23% | 2.04757003794e-19 |
| 88 | 1.75921860444e+13 · [Bos–Conway–Sloane BCS88 (1982)](https://neilsloane.com/doc/Me103.pdf) | 6.58676184812e+13 · [Primitive MW92 kernel section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | +274.41% | 1.85638056287e-19 |
| 89 | 7.14465044353e+13 · [Guaranteed MW92 section via Rogers’ inequality](research/covering_and_deep_cosets_20261115/SECTIONS_PROOF.md) (derived) | 2.03871499742e+14 · [Primitive MW92 kernel section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | +185.27% | 1.52239236465e-19 |
| 90 | 3.61283064288e+14 · [Guaranteed MW92 section via Rogers’ inequality](research/covering_and_deep_cosets_20261115/SECTIONS_PROOF.md) (derived) | 7.28635694099e+14 · [Primitive MW92 kernel section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | +101.64% | 1.43364905437e-19 |
| 91 | 1.86165981730e+15 · [Guaranteed MW92 section via Rogers’ inequality](research/covering_and_deep_cosets_20261115/SECTIONS_PROOF.md) (derived) | 2.48346380530e+15 · [Primitive MW92 kernel section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | +33.38% | 1.28046110298e-19 |
| 92 | 9.77405088148e+15 · [Shioda MW92](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | — | — | 1.31340615344e-19 |
| 93 | 4.88702544074e+15 · [Shioda MW92 x classical 1D lattice](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) (derived) | 8.74303398226e+15 · [MW92 lattice shear](research/mw92_extensions_20261122/README.md#lattice-shears) (existence) | +78.90% | 3.04556337848e-20 |
| 94 | 2.82152545374e+15 · [Shioda MW92 x classical 2D lattice](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) (derived) | 8.94634610958e+15 · [MW92 order-three dichotomy](research/mw92_order3_dichotomy_20261122/README.md) (existence) | +217.07% | 8.03567126556e-21 |
| 95 | 1.72782441448e+15 · [Shioda MW92 x classical 3D lattice](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) (derived) | 9.13261126748e+15 · [MW92 order-three dichotomy](research/mw92_order3_dichotomy_20261122/README.md) (existence) | +428.56% | 2.10405317054e-21 |
| 96 | 4.75445050459e+15 · [BCS96, Table 1 (1982)](https://neilsloane.com/doc/Me103.pdf) | 9.60319359100e+15 · [Refined MW92 order-three lattice](research/mw92_96_equal_cover_20261122/README.md) (existence) | +101.98% | 5.64547885516e-22 |

## Scope and replay

- The 55D thesis target is disputed: the displayed lattice has an explicit shorter vector. The row compares with the [rigorous repair](research/tilted55_20261116/PROOF.md), not the printed target. It is excluded from the improvement count.
- The 89–91D comparison intervals are bounds for a cited derived construction guarantee, not upper bounds on the optimal density. The same distinction applies to historical printed decimal prefixes in other rows.
- The 79D and 86–91D periodic packings and the nine-translate 96D packing are not counted as lattices. In 86–91D the listed bare kernels form single lattices and retain the parent minimum lower bound.
- This catalog makes no lattice-optimality claim. Printed and derived comparisons, exact squared values, and the proof paths are in [catalog.json](research/lattice_lower_49_96_20261122/catalog.json).

Verified improvements over the named non-disputed prior bounds: **31 dimensions**.

Run the Python script in [the proof directory](research/lattice_lower_49_96_20261122/verify_catalog.py) to check the table, or add --write to rebuild it.
