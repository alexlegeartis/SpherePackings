# Sphere-packing density bounds and constructions

[Henry Cohn’s sphere-packing tables](https://cohn.mit.edu/sphere-packing/) cite this repository for the following records (**checked 17 November 2026**):
- **Unrestricted upper bounds:** dimensions **4–7, 9–23 and 25–48**.
- **Lattice-only upper bounds:** dimensions **14–23 and 25–48**.
- **Packing lower bound:** dimension **36**, in both tables.

[Verified citation coverage and source snapshot](research/cohn_listing_20261117/). The current page lists dimensions 1–48; this is not a claim that every latest refinement, or any 49D+ result, is already listed. Citation does not establish optimality or independent verification by the table’s author.

Reproducible certificates for unrestricted sphere packings, lattice-only upper bounds and packing constructions are separated below. Comparisons remain source-attributed; derived guarantees, numerical benchmarks and certified bounds are not interchangeable, and no exhaustive worldwide-priority claim is made.

[General lower bounds](#general-lower-bounds) · [General upper bounds](#general-upper-bounds) · [Lattice-only lower bounds](#lattice-only-lower-bounds) · [Lattice-only upper bounds](#lattice-only-upper-bounds)

<a id="packing-constructions-lower-bounds"></a>

## General lower bounds

[All occupied-density lower bounds in dimensions 1–96](LOWER_BOUNDS_1_96.md).

**23 November:** [An exact 96D four-Leech family ceiling](research/gamma96_fourleech_motif96_20261123/) rules out beating the current periodic bound when the three-block Γ72 carrier and integral horizontal phases are retained and the fourth Leech block has arbitrary scalar scale. The proof covers arbitrary nonlinear motifs and multiple fibres at a vertical residue. It does not exclude larger phase alphabets or arbitrary four-block geometry.

**23 November:** [A quantitative 39D obstruction](research/e7_half_radical39_20261123/) rules out improving the current bound by assigning the five-E7 half-radical phases to every site of any 4D vertical lattice at the inherited minimum. It allows nonlinear phase assignment and small perturbations of phases and vertical sites. In 48D, [ternary phases cannot attain one-period equality](research/leech_order3_equality48_20261123/) and [mixed half-period motifs with a second factor of minimum at most 3](research/leech_mixed_half_period48_20261123/) stay below 72% of P48. These are family restrictions, not global upper bounds.

**23 November:** [Quantitative metric rigidity](research/p48_gamma72_metric_and_equality_20261123/) proves that no arbitrary metric of condition number at most 2 improves a fixed extremal 48D even unimodular lattice. The corresponding 72D threshold is 19/10. The same package classifies equality in the 48D additive Leech half-period family by a binary quadratic-form isometry. These are scoped statements, not global density bounds.

**23 November:** [Two-parity D3 stacking](research/e7_periodic38_20261123/) reuses the E7 glue to raise the **38D unrestricted lower bound by 19.33%**. It is a nonlattice periodic packing. A separate [coupled E7 lattice lift](research/e7_lattice38_20261123/) now attains the same density and improves the preceding **lattice-only** bound by **19.33%**. Its [primitive 37D section](research/e7_37_section_20261123/) improves the preceding lattice and unrestricted 37D bound by a further **15.47%**.

**23 November:** [One five-block E7 construction](research/e7_hermitian35_20261123/) raises the lattice and unrestricted density lower bounds in **34–37D** by **3.93%, 41.42%, 51.47%, and 41.42%**, respectively. The 36D comparison is the preceding repository result listed by Cohn. The proof also gives an exact optimum within the fixed 35D additive binary-phase family.

**23 November:** [Degree-13 covering certificates](research/mw92_cfr_taylor13_20261123/) sharpen the 93–96D unrestricted density bounds, including the lattice bounds in 93–95D. Gains over the preceding certificates are all below 0.005%; the 96D lattice bound is unchanged.

**23 November:** [Exact odd-Leech product-period optimum](research/odd_leech_product48_20261123/) matches P48 for arbitrary half-period motifs with two orthogonal, independently scaled odd-Leech factors.

**23 November:** [Ternary phase-depth bounds](research/order_three_tori_48_72_20261123/) cap the squared depth of every order-three torus phase by 4 in any extremal 48D even unimodular lattice and 6 in any extremal 72D one. The proof also gives sharper residue-dependent bounds.

**23 November:** [One-Leech-period projection bounds](research/leech_single_period_20261123/) give exact scoped optima in 48D and in 25–36D for lattices (symbolic in 34–36D), including arbitrary periodic motifs in 25–27D and 32D. The result allows arbitrary transverse periods and applies in the stated high-minimum regimes.

**23 November:** [Exact 48D scalar-shear bound](research/leech_scalar_shear48_20261123/) extends the Leech product-period optimum to arbitrary real scalar coupling between the two factors. It covers every half-period motif and arbitrary real shifts at normalized squared minimum at least 6. It does not settle general 48D packing.

**23 November:** [Leech product-period optima](research/leech_product_periods_48_72_20261123/) match the 48D and 72D bounds in scoped families of periodic packings. The 48D result includes arbitrary real shifts when the squared minimum is at least 6 in the stated normalization. These are not global upper bounds.

**23 November:** [A six-translate Ne64 packing](research/ne64_dual_matching59_20261123/) improves the unrestricted 59D density lower bound by **0.3246%**. It is not a lattice.

**23 November:** [A four-layer shear](research/vardy30_shear31_20261123/) raises the unrestricted 31D packing lower bound by **16.92%**. The lattice-only bound is unchanged.

**22 November:** [Half-root Γ72–D4 lattices](research/gamma72_root_tilt_20261122/) improve the preceding 77D and 78D bounds by **8.01%** and **21.72%**, respectively (net **15.47%** and **33.33%** over the bounds before the tilt search). The historical 79D lattice has since been superseded by the [bare MW80 section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#79d-bare-mw80-section). The denser 79D periodic packing remains the unrestricted bound.

**22 November:** six explicit lattices improve the preceding 73–78D lattice bounds by **11.80% to 56.25%**. The 78D lattice also improves the preceding unrestricted bound by **34.59%**. [Exact construction and comparisons](research/gamma72_depth6_20261122/).

**17 November binary-CM update:** certified 57D **≥1.61599884109140e-11** and 58D **≥5.58243739609337e-12**. [Exact parent, complete independent coset exclusion and extension proof](research/packing_binary57_58_20261117/). Comparisons retain the pre-existing published-source consequences, not our own preceding constructions.

**17 November update: six independently audited periodic packing improvements in 85D–90D.** The table compares with attributable external constructions and certified consequences of published theorems, not our own preceding results. [Exact construction and replay](research/packing_release8590_20261117/) · [External attribution, exact gain checks and separate internal progress](research/external_comparison85_90_20261117/). Previous constructions remain available.

**Comparison method / source** names the published work behind each recorded
comparison. Values marked **derived** are calculated from published inputs,
rather than quoted as numerical records in those sources. All densities in
this table are occupied-volume fractions.

The [literature review dated 16 November 2026](research/covering_and_deep_cosets_20261115/PRIORITY.md) records the external guarantees and derived comparisons. Attribution and exact comparison arithmetic for 85D–90D were rechecked on 17 November. **Derived** does not mean that the source printed that dimension-specific record. For 89D–90D, the comparison is a newly evaluated consequence of pre-existing spherical linear-programming and section theorems, not a numerical record printed by those authors. It improves on the earlier Rogers-derived comparison. These are the strongest external-derived comparisons certified in this follow-up, not an exhaustive worldwide-priority claim. Numerical bootstrap spectra suggest stronger expected Cohn–Elkies (CE) comparisons for the 89D–90D section route, but the same-radius lift to a CE function is conjectural in general and the required function inequalities are not certified here; see the linked comparison note.

| n | Comparison density / method | Our density / method | Increase |
| ---: | --- | --- | ---: |
| 29 | ≈3.41446469074224e-5<br>[B29* — Vardy (1999)](https://doi.org/10.1090/S0002-9947-99-02169-8) · [Cohn’s unrestricted table](https://cohn.mit.edu/sphere-packing/#data) | **≥4.82878227389174e-5**<br>[Alternating Vardy-type layers, exact center density `1`](research/vardy28_stack29_20261122/) | **41.421356%** |
| 31 | ≈1.18377651859338e-5<br>[Q31 — Quebbemann (1984), section of Q32](https://doi.org/10.1112/S0025579300010731) · [Cohn’s unrestricted table](https://cohn.mit.edu/sphere-packing/#data) | **≥1.38411060495194e-5**<br>[Order-four shear of Vardy-type stocks, exact center density `sqrt(2)`](research/vardy30_shear31_20261123/) | **16.923303%** |
| 34 | ≈1.76697388915406e-6<br>[Cohn’s pinned table, ref. 13](https://cohn.mit.edu/sphere-packing/#13) | **≥1.83644820734493e-6**<br>[Hermitian E7 glue: primitive section; center density `4/sqrt(3)`](research/e7_hermitian35_20261123/) | **3.93%** |
| 35 | ≈9.46190415115372e-7<br>[Cohn’s pinned table, ref. 13](https://cohn.mit.edu/sphere-packing/#13) | **≥1.33811531764358e-6**<br>[Five-block Hermitian E7 glue; center density `4`](research/e7_hermitian35_20261123/) | **41.42%** |
| 36 | ≈7.33034778316707e-7<br>[Previous 36D metric-deformation lattice](constructions/dim36/) (listed by Cohn) | **≥1.11031619697710e-6**<br>[Hermitian E7 glue: order-four lift; center density `8`](research/e7_hermitian35_20261123/) | **51.47%** |
| 37 | ≈3.21356200759300e-7<br>[Cohn’s pinned table, ref. 13](https://cohn.mit.edu/sphere-packing/#13) | **≥5.24772478359786e-7**<br>[Primitive 37D section of the coupled E7 lift; center density `16/sqrt(3)`](research/e7_37_section_20261123/) | **63.30%** |
| 38 | ≈2.17566730891289e-7<br>[Previous antipode construction (2026)](https://arxiv.org/abs/2607.20359) · [Cohn’s unrestricted table](https://cohn.mit.edu/sphere-packing/#data) | **≥2.59631836184760e-7**<br>[Coupled E7 lattice, exact center density `8*sqrt(2)`](research/e7_lattice38_20261123/); [earlier periodic construction](research/e7_periodic38_20261123/) | **19.33%** |
| 49 | ≈4.128837170599e-9 (**derived**)<br>P48q × A1 — [Bos–Conway–Sloane (1982), §5](https://neilsloane.com/doc/Me103.pdf); [exact comparison](research/covering_and_deep_cosets_20261115/literature_comparisons.json) | ≈7.151355755657e-9<br>[Explicit deeper P48q half-coset; full minimum certificate](research/covering_and_deep_cosets_20261115/DEEPER_PROOF.md) | **73.205081%** |
| 50 | ≈8.408155459659e-10 (**derived**)<br>P48q × A2 — [Bos–Conway–Sloane (1982), §5](https://neilsloane.com/doc/Me103.pdf); [exact comparison](research/covering_and_deep_cosets_20261115/literature_comparisons.json) | ≈2.522446637898e-9<br>[Deep plane + fitted metric; explicit minimum-6 lattice](research/adaptive_layers_20261115/P48_PROOF.md) | **200.000000%** |
| 51 | ≈1.798425511734e-10 (**derived**)<br>P48q × A3 — [Bos–Conway–Sloane (1982), §5](https://neilsloane.com/doc/Me103.pdf); [exact comparison](research/covering_and_deep_cosets_20261115/literature_comparisons.json) | **≥9.34489307985541e-10**<br>[Explicit P48n quaternionic lattice](research/p48n_hurwitz_20261122/) | **419.615242%** |
| 52 | ≈1.094979637686e-10 (**derived**)<br>Two Blichfeldt sections of Leterrier54 — [Leterrier (2022), Corollary 1.4](https://doi.org/10.1007/s40993-022-00321-0); [derivation](research/range48_96_20261115/verify_range.py) | **≥3.95931817049999e-10**<br>[Explicit P48n quaternionic lattice](research/p48n_hurwitz_20261122/) | **261.588292%** |
| 53 | ≈1.256125809644e-10 (**derived**)<br>One Blichfeldt section of Leterrier54 — [Leterrier (2022), Corollary 1.4](https://doi.org/10.1007/s40993-022-00321-0); [comparison derivation](research/range48_96_20261115/verify_range.py) | **≥1.66176494990127e-10**<br>[Explicit P48n quaternionic lattice with a fifth layer](research/p48n_hurwitz_20261122/) | **32.292875%** |
| 57 | ≈8.85119018168e-12 (**derived**)<br>MW56 × A1 — [Schütt–Shioda, Theorem 13.35, p=29](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | **≥1.61599884109140e-11**<br>[Exact depth-7 Mordell–Weil half-coset extension](research/packing_binary57_58_20261117/) | **82.574186%** |
| 58 | ≈4.474935526505e-12 (**derived**)<br>Γ72 perpendicular to a rank-14 Leech block — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [exact derived section](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | **≥5.58243739609337e-12**<br>[Exact depth-7 Mordell–Weil half-coset extension with complex-multiplication triangular layers](research/packing_binary57_58_20261117/) | **24.749002%** |
| 59 | ≈2.585163653787e-12 (**derived**)<br>Γ72 perpendicular to √2κ13 — [Nebe (2012), §4.2](https://arxiv.org/pdf/1008.2862) + [classical κ13 catalogue data](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KAPPA13.html); [derivation](constructions/dim59/#literature-derived-comparison) | **≥2.64616393271302e-12**<br>[Six-translate Ne64 antipode with matched normal](research/ne64_dual_matching59_20261123/) | **2.359629%** |
| 60 | ≈[1.668964977674e-12, 1.680573547119e-12)<br>Elkies section guarantee — [Vardy (1999), §4, p.278](https://doi.org/10.1090/S0002-9947-99-02169-8); published log₂ δ prefix 19.04… | ≈2.395250641682e-12<br>[Explicit Ne64 antipode; exact period and shifts](research/novel_routes_20261113/audits/README.md) | **>42.525785%** |
| 61 | ≈[1.508819950024e-12, 1.519314622713e-12)<br>Elkies section guarantee — [Vardy (1999), §4, p.278](https://doi.org/10.1090/S0002-9947-99-02169-8); published log₂ δ prefix 20.54… | ≈2.022035232382e-12<br>[Explicit Ne64 antipode; exact period and shifts](research/novel_routes_20261113/audits/README.md) | **>33.088644%** |
| 62 | ≈[1.440184355054e-12, 1.450201629426e-12)<br>Elkies section guarantee — [Vardy (1999), §4, p.278](https://doi.org/10.1090/S0002-9947-99-02169-8); published log₂ δ prefix 22.13… | ≈1.628191894214e-12<br>[Explicit Ne64 antipode; exact period and shifts](research/novel_routes_20261113/audits/README.md) | **>12.273483%** |
| 65 | ≈3.201378537714e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 7 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈3.896264824634e-13<br>[Explicit Ne64 extension; complete minimum-12 certificate](research/novel_routes_20261113/higher_parent/) | **21.705846%** |
| 66 | ≈2.272528652057e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 6 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈2.431106576448e-13<br>[Explicit Γ72 antipode; exact normal, period and shifts](research/novel_routes_20261113/higher_parent/) | **6.978039%** |
| 67 | ≈1.698310457767e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 5 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈2.215613494306e-13<br>[Explicit Γ72 antipode; exact normal, period and shifts](research/novel_routes_20261113/higher_parent/) | **30.459863%** |
| 68 | ≈1.454792126225e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 4 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈2.200710220661e-13<br>[Explicit Γ72 antipode; exact normal, period and shifts](research/novel_routes_20261113/higher_parent/) | **51.273174%** |
| 69 | ≈1.23719412351e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 3 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈1.908606667819e-13<br>[Explicit Γ72 antipode; exact normal, period and shifts](research/novel_routes_20261113/higher_parent/) | **54.268973%** |
| 70 | ≈1.206263329664e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 2 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈1.73198888183e-13<br>[Explicit Γ72 antipode; exact normal, period and shifts](research/novel_routes_20261113/higher_parent/) | **43.582984%** |
| 71 | ≈1.238694909866e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 1 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈1.417772488475e-13<br>[Explicit Γ72 antipode; exact normal, period and shifts](research/novel_routes_20261113/higher_parent/) | **14.456956%** |
| 73 | ≈2.132314974725e-14 (**derived**)<br>Γ72 × A1 — [Nebe (2012)](https://arxiv.org/abs/1008.2862) | **≥4.26462994945094e-14**<br>[Exact Γ72 half-coset lattice](research/gamma72_depth6_20261122/) | **100%** |
| 74 | ≈3.575177737935e-15 (**derived**)<br>Γ72 × A2 — [Nebe (2012)](https://arxiv.org/abs/1008.2862) | **≥1.43007109517410e-14**<br>[Exact Γ72 half-coset lattice](research/gamma72_depth6_20261122/) | **300%** |
| 75 | ≈6.315750546658e-16 (**derived**)<br>Γ72 × A3 — [Nebe (2012)](https://arxiv.org/abs/1008.2862) | **≥5.05260043732660e-15**<br>[Exact Γ72 half-coset lattice](research/gamma72_depth6_20261122/) | **700%** |
| 76 | ≈1.279866139755e-16 (**derived**)<br>Γ72 × D4 — [Nebe (2012)](https://arxiv.org/abs/1008.2862) | **≥2.04778582360778e-15**<br>[Exact Γ72 half-coset lattice](research/gamma72_depth6_20261122/) | **1500%** |
| 77 | ≈9.521458141249e-17 (**derived**)<br>Bare L80 section, normal rank 3 — [Bachoc–Nebe (1998)](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf); [normal matrix audit](research/cross_dimension_lower_20261115/l80_results.json) | **≥6.73268761846130e-16**<br>[Half-root Γ72–D4 lattice](research/gamma72_root_tilt_20261122/) | **607.106781%** |
| 78 | ≈9.523775961404e-17 (**pinned external-derived example**)<br>L80 perpendicular to2A2 — [Bachoc–Nebe(1998)](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf) · [Exact derivation and scope](research/packing78_verified_20261117/) | **≥2.53967358970778e-16**<br>[Half-root Γ72–D4 lattice](research/gamma72_root_tilt_20261122/) | **166.666666%** |
| 79 | ≈9.274742674901e-17 (**derived**)<br>L80 perpendicular to a minimum vector — [Bachoc–Nebe (1998)](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf) | ≈1.119208813283e-16<br>[Explicit MW80 antipode; exact period, shifts and minimum attainment](research/creative_sections_20261116/PROOF.md) | **20.672762%** |
| 81 | ≈1.595311593607e-17 (**derived**)<br>MW80 × A1 — [Schütt–Shioda (2017 author draft), Theorem 13.35, p=41](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf); [exact comparison](research/covering_and_deep_cosets_20261115/literature_comparisons.json) | **≥2.98455470421258e-17**<br>[MW80 half-coset lift](research/mw80_half_coset_81_20261123/README.md) | **87.082869%** |
| 82 | ≈2.541814578499e-18 (**derived**)<br>MW80 × A2 — [Schütt–Shioda, Theorem 13.35, p=41](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | **≥9.17033181889404e-18**<br>[All-deep A2 section](research/l80_a2_section_20261122/README.md) | **260.778944%** |
| 83 | ≈5.446261557232e-19 (**derived**)<br>BCS88 perpendicular to √2 D5 — [Bos–Conway–Sloane (1982), Theorems 1–2, Eq. (15)](https://neilsloane.com/doc/Me103.pdf); [exact section audit](research/bcs83_comparator_20261117/) | **≥3.08087078338732e-18**<br>[Bachoc–Nebe lattice + independently audited overlapping deep-coset parity extension](research/packing83_85_deep_verified_20261117/README_DEEP.md) | **>465.685424%** |
| 84 | ≈2.970203932666e-19 (**derived**)<br>Section of BCS88 perpendicular to √2 D4 — [Bos–Conway–Sloane (1982), Theorems 1–2, Eq. (15)](https://neilsloane.com/doc/Me103.pdf); [derivation and exact matrix audit](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) | **≥1.18808157306620e-18**<br>[Bachoc–Nebe lattice + independently audited overlapping deep-coset parity extension](research/packing83_85_deep_verified_20261117/README_DEEP.md) | **300%** |
| 85 | ≈1.610346851551e-19 (**derived**)<br>Section of BCS88 perpendicular to √2 A3 — [Bos–Conway–Sloane (1982), Theorems 1–2, Eq. (15)](https://neilsloane.com/doc/Me103.pdf); [derivation and exact matrix audit](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) | **≥4.55474871517601e-19**<br>[Bachoc–Nebe lattice + independently audited overlapping deep-coset parity extension](research/packing83_85_deep_verified_20261117/README_DEEP.md) | **>182.842712%** |
| 86 | ≈1.002298419482e-19 (**derived**)<br>BCS88 perpendicular to √2A2 — [Bos–Conway–Sloane (1982), Theorems 1–2](https://neilsloane.com/doc/Me103.pdf); [derived section](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) | **≥3.2157847626326119e-19**<br>[MW92 projected-simplex packing](research/mw92_small_simplexes_20261122/README.md) | **>220.8410%** |
| 87 | ≈6.578926483929e-20 (**derived**)<br>BCS88 perpendicular to √2A1 — [Bos–Conway–Sloane (1982), Theorems 1–2](https://neilsloane.com/doc/Me103.pdf); [derived section](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) | **≥2.6231014000959713e-19**<br>[MW92 projected-simplex packing](research/mw92_small_simplexes_20261122/README.md) | **>298.7126%** |
| 88 | ≈4.95809518917e-20<br>Minimum-8 unimodular lattices — [Sakamoto (2026)](https://arxiv.org/abs/2609.15633v2); same density as [BCS1982](https://neilsloane.com/doc/Me103.pdf) | **≥2.42250297455694e-19**<br>[Correlated Mordell–Weil dual-normal sections with optimized simplex translate clusters](research/packing_release8590_20261117/PROOF.md) | **>388.59%** |
| 89 | ≈5.435520207304e-20 (**derived**)<br>MW92 sections via Levenshtein spherical LP and [Cohn–Zhao (2014), Proposition2.1](https://arxiv.org/html/1212.5966) — parent: [Schütt–Shioda, Theorem13.35](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf); [exact derivation and independent check](research/external_comparison85_90_20261117/KL_PROOF.md) | **≥1.98411302345717e-19**<br>[Correlated Mordell–Weil dual-normal sections with optimized simplex translate clusters](research/packing_release8590_20261117/PROOF.md) | **>265.02%** |
| 90 | ≈7.200357004613e-20 (**derived**)<br>MW92 sections via Levenshtein spherical LP and [Cohn–Zhao (2014), Proposition2.1](https://arxiv.org/html/1212.5966) — parent: [Schütt–Shioda, Theorem13.35](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf); [exact derivation and independent check](research/external_comparison85_90_20261117/KL_PROOF.md) | **≥1.61584267273781e-19**<br>[Correlated Mordell–Weil dual-normal sections with optimized simplex translate clusters](research/packing_release8590_20261117/PROOF.md) | **>124.41%** |
| 91 | ≈[9.598621803784e-20, 9.59946031323e-20] (**derived**)<br>MW92 / full Rogers-simplex section — [Schütt–Shioda, Theorem 13.35](<https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf>) + [Rogers (1958)](https://doi.org/10.1112/plms/s3-8.4.609); [derived interval](research/mw92_sections_20261115/PROOF.md#external-comparisons-including-the-full-rogers-consequence) | **≥1.28046110298280e-19**<br>[Primitive MW92 bare lattice section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | **>33.38%** |
| 93 | ≈1.702354782363e-20 (**derived**)<br>MW92 × A1 — [Schütt–Shioda, Theorem 13.35, p=47](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | **≥3.045564894259504e-20**<br>[MW92 degree-13 covering refinement; lattice existence proof](research/mw92_cfr_taylor13_20261123/) | **>78.9030%** |
| 94 | ≈2.534314091589e-21 (**derived**)<br>MW92 × A2 — [Schütt–Shioda, Theorem 13.35, p=47](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | **≥8.036027386567117e-21**<br>[MW92 degree-13 covering refinement; lattice existence proof](research/mw92_cfr_taylor13_20261123/) | **>217.0888%** |
| 95 | ≈3.980717377502e-22 (**derived**)<br>MW92 × A3 — [Schütt–Shioda, Theorem 13.35, p=47](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | **≥2.104056220143222e-21**<br>[MW92 degree-13 covering refinement; lattice existence proof](research/mw92_cfr_taylor13_20261123/) | **>428.5620%** |
| 96 | ≈2.795023294827e-22<br>P48q-based lattice — [Bos–Conway–Sloane (1982), Eq. (15), Table 1, p.177](https://neilsloane.com/doc/Me103.pdf); exact center density `3^48/2^24` | **≥6.250833919802144e-22**<br>[MW92 degree-13 covering refinement; nine-translate periodic existence proof](research/mw92_cfr_taylor13_20261123/) | **>123.6415%** |

**Scope:** Some rows are existence bounds, not explicit packings. Named comparisons include derived guarantees; they do not establish exhaustive priority. [Sources and exact gains](research/covering_and_deep_cosets_20261115/PRIORITY.md) · [Full audit and research notes](RELEASE_NOTES_20261117.md#packing-comparison-and-research-notes).

**55D** now has a verified occupied density of **≥7.54599036285564e-11**, a **0.9452%** increase over the preceding repository bound. The center density is `15/19` of the full printed challenge target. [Three-layer construction and exact check](research/uneven_layers55_20261122/README.md) · [Printed comparison audit](research/cross_coset55_20261116/PROOF.md#external-comparison).

<a id="registered-results"></a>

## General upper bounds

**Recorded upper U** and **lower L** are the registered reference values used for the percentages in that row. For current lower bounds, use the [lower-bound catalog](LOWER_BOUNDS_1_96.md). **Comparison method / source** names the method behind U and links to its paper or pinned certificate; **Our method / proof** describes the new bound N. A recorded comparison is not a claim about the latest worldwide record. [Full source titles, versions and comparison scope](COMPARISONS.md#table-comparison-methods-and-sources).

| Dimension / lower L | Recorded upper U / source | Our upper N / method | Upper decrease | Gap decrease |
|---|---|---|---:|---:|
| 4<br>**L:** [≈0.61685](https://cohn.mit.edu/sphere-packing/) | 0.629521681653<br>[Scalar/kernel + count and far-pair refinements — Felpix-Studios, pinned bb97550](https://github.com/Felpix-Studios/sphere-packing/tree/bb97550ba47ad92bacfee567843888be16895a23/preprints/D4) | <0.629520189662<br>[Degree16 far-pair polynomial + unused angular Gram credits](verifications/dim4/refinement_energy20261115/) · [Sources & formulas](COMPARISONS.md) | 0.000237% | 0.011774% |
| 5<br>**L:** [≈0.46526](https://cohn.mit.edu/sphere-packing/) | ≈0.51226<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.50995<br>[Spherical lens caps + long-edge PSD accounting](verifications/dim5/refinement_spherical_lens/) · [Sources & formulas](COMPARISONS.md) | 0.45% | 4.92% |
| 6<br>**L:** [≈0.37295](https://cohn.mit.edu/sphere-packing/) | ≈0.41029<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.40913<br>[Neighbor-count cuts + Schur-saturated cyclic kernel](verifications/dim6/refinement_patch5/) · [Sources & formulas](COMPARISONS.md) | 0.28% | 3.11% |
| 7<br>**L:** [≈0.29530](https://cohn.mit.edu/sphere-packing/) | ≈0.32111<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.32058<br>[Neighbor-count cuts + Schur-saturated cyclic kernel](verifications/dim7/) · [Sources & formulas](COMPARISONS.md) | 0.17% | 2.05% |
| 9<br>**L:** [≈0.14577](https://cohn.mit.edu/sphere-packing/) | ≈0.19099<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.19093<br>[Signed-hinge cyclic profile + average-neighbor count](verifications/dim9/) · [Sources & formulas](COMPARISONS.md) | 0.03% | 0.14% |
| 10<br>**L:** [≈0.099616](https://cohn.mit.edu/sphere-packing/) | ≈0.14340<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.14333<br>[Log-hinge cyclic profile + root-cut radial repair](verifications/dim10/) · [Sources & formulas](COMPARISONS.md) | 0.05% | 0.17% |
| 11<br>**L:** [≈0.066238](https://cohn.mit.edu/sphere-packing/) | ≈0.10672<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.10664<br>[Quartic radial damping + all-dilation cyclic reduction](verifications/dim11/) · [Sources & formulas](COMPARISONS.md) | 0.08% | 0.21% |
| 12<br>**L:** [≈0.049454](https://cohn.mit.edu/sphere-packing/) | ≈0.079707<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.079503<br>[Quartic radial damping + all-dilation cyclic reduction](verifications/dim12/) · [Sources & formulas](COMPARISONS.md) | 0.26% | 0.67% |
| 13<br>**L:** [≈0.032014](https://cohn.mit.edu/sphere-packing/) | ≈0.060121<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.059905<br>[Quartic radial damping + all-dilation cyclic reduction](verifications/dim13/) · [Sources & formulas](COMPARISONS.md) | 0.36% | 0.77% |
| 14<br>**L:** [≈0.021624](https://cohn.mit.edu/sphere-packing/) | ≈0.045010<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.045003<br>[Quartic radial damping + all-dilation cyclic reduction](verifications/dim14/) · [Sources & formulas](COMPARISONS.md) | 0.02% | 0.03% |
| 15<br>**L:** [≈0.016858](https://cohn.mit.edu/sphere-packing/) | ≈0.033726<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.033672<br>[Quartic radial damping + all-dilation cyclic reduction](verifications/dim15/) · [Sources & formulas](COMPARISONS.md) | 0.16% | 0.32% |
| 16<br>**L:** [≈0.014708](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/BW16.html) | ≈0.024921<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.024901<br>[Widened shell count + coupled log-hinge radial kernel](verifications/dim16/widened_reconstruction/) · [Sources & formulas](COMPARISONS.md) | 0.08% | 0.20% |
| 17<br>**L:** [≈0.0088113](https://doi.org/10.4153/CJM-1964-065-1) | ≈0.018464<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.018035<br>[Broader Gaussian radial profile + degree95 Fourier enrichment](verifications/dim17/refinement_broad20261118/) · [Sources & formulas](COMPARISONS.md) | 2.32% | 4.44% |
| 18<br>**L:** [≈0.0061679](https://doi.org/10.1023/A:1008733715001) | ≈0.013485<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.012945<br>[Broader Gaussian radial profile + degree95 Fourier enrichment](verifications/dim18/refinement_broad20261118/) · [Sources & formulas](COMPARISONS.md) | 4.01% | 7.39% |
| 19<br>**L:** [≈0.0041474](https://arxiv.org/abs/2505.02394v1) | ≈0.0098180<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0092288<br>[Broader Gaussian radial profile + degree95 Fourier enrichment](verifications/dim19/refinement_broad20261118/) · [Sources & formulas](COMPARISONS.md) | 6.00% | 10.39% |
| 20<br>**L:** [≈0.0040242](https://arxiv.org/abs/2505.02394v1) | ≈0.0071271<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0066416<br>[Broader Gaussian radial profile + degree95 Fourier enrichment](verifications/dim20/refinement_broad20261118/) · [Sources & formulas](COMPARISONS.md) | 6.81% | 15.65% |
| 21<br>**L:** [≈0.0029298](https://arxiv.org/abs/2505.02394v1) | ≈0.0051597<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0048577<br>[Broader Gaussian radial profile + degree95 Fourier enrichment](verifications/dim21/refinement_broad20261118/) · [Sources & formulas](COMPARISONS.md) | 5.85% | 13.54% |
| 22<br>**L:** [≈0.0024510](https://cohn.mit.edu/sphere-packing/) | ≈0.0037259<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0035915<br>[Broader Gaussian radial profile + degree95 Fourier enrichment](verifications/dim22/refinement_broad20261118/) · [Sources & formulas](COMPARISONS.md) | 3.61% | 10.54% |
| 23<br>**L:** [≈0.0019072](https://cohn.mit.edu/sphere-packing/) | ≈0.0026843<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0026052<br>[Broader Gaussian radial profile + degree95 Fourier enrichment](verifications/dim23/refinement_broad20261118/) · [Sources & formulas](COMPARISONS.md) | 2.95% | 10.18% |
| 25<br>**L:** [≈0.00067721](https://cohn.mit.edu/sphere-packing/) | ≈0.0013842<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0013611<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim25/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 1.67% | 3.27% |
| 26<br>**L:** [≈0.00026922](https://cohn.mit.edu/sphere-packing/) | ≈0.00099102<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00096119<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim26/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 3.01% | 4.13% |
| 27<br>**L:** [≈0.00015759](https://cohn.mit.edu/sphere-packing/) | ≈0.00070823<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00067374<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim27/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 4.87% | 6.26% |
| 28<br>**L:** [≈0.00010464](https://cohn.mit.edu/sphere-packing/) | ≈0.00050525<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00047150<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim28/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 6.68% | 8.43% |
| 29<br>**L:** [≈0.000048288](research/vardy28_stack29_20261122/) | ≈0.00035986<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00033305<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim29/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 7.45% | 8.60% |
| 30<br>**L:** [≈0.000021915](https://cohn.mit.edu/sphere-packing/) | ≈0.00025590<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00023903<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim30/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 6.59% | 7.21% |
| 31<br>**L:** [≈0.000011838](https://cohn.mit.edu/sphere-packing/) | ≈0.00018171<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00017253<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim31/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 5.05% | 5.40% |
| 32<br>**L:** [≈0.000011041](https://cohn.mit.edu/sphere-packing/) | ≈0.00012884<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00012276<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim32/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 4.72% | 5.16% |
| 33<br>**L:** [≈0.0000041407](https://cohn.mit.edu/sphere-packing/) | ≈0.000091236<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000089087<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim33/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 2.35% | 2.47% |
| 34<br>**L:** [≈0.0000017670](https://cohn.mit.edu/sphere-packing/) | ≈0.000064522<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000063254<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim34/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 1.97% | 2.02% |
| 35<br>**L:** [≈9.4619e-7](https://cohn.mit.edu/sphere-packing/) | ≈0.000045574<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000043911<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim35/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 3.65% | 3.73% |
| 36<br>**L:** [≈7.3303e-7 (this work)](constructions/dim36/) | ≈0.000032153<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000030495<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim36/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 5.16% | 5.28% |
| 37<br>**L:** [≈3.2136e-7](https://cohn.mit.edu/sphere-packing/) | ≈0.000022659<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000020943<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim37/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 7.57% | 7.68% |
| 38<br>**L:** [≈2.1757e-7](https://arxiv.org/html/2607.20359v4) | ≈0.000015951<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000014518<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim38/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 8.98% | 9.10% |
| 39<br>**L:** [≈1.6803e-7](https://arxiv.org/html/2607.20359v4) | ≈0.000011217<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000010220<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim39/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 8.89% | 9.02% |
| 40<br>**L:** [≈1.4798e-7](https://arxiv.org/html/2607.20359v4) | ≈0.0000078801<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0000073057<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim40/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 7.29% | 7.43% |
| 41<br>**L:** [≈8.6368e-8](https://arxiv.org/html/2607.20359v4) | ≈0.0000055306<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0000051253<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim41/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 7.33% | 7.45% |
| 42<br>**L:** [≈5.7517e-8](https://arxiv.org/html/2607.20359v4) | ≈0.0000038781<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0000036196<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim42/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 6.67% | 6.77% |
| 43<br>**L:** [≈4.6423e-8](https://arxiv.org/html/2607.20359v4) | ≈0.0000027169<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0000026096<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim43/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 3.95% | 4.02% |
| 44<br>**L:** [≈3.9244e-8](https://arxiv.org/abs/2505.02394v1) | ≈0.0000019018<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0000018326<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim44/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 3.64% | 3.71% |
| 45<br>**L:** [≈3.5582e-8](https://arxiv.org/abs/2505.02394v1) | ≈0.0000013301<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0000012613<br>[Reoptimized nonnegative Gegenbauer cyclic kernel + signed logarithmic profile + retained Fourier witness](verifications/general_angular_20261117/) · [Sources & formulas](COMPARISONS.md) | 5.17% | 5.31% |
| 46<br>**L:** [≈2.8610e-8](https://doi.org/10.1007/BF01232378) | ≈9.2952e-7<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈8.7124e-7<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim46/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 6.27% | 6.47% |
| 47<br>**L:** [≈2.2670e-8](https://arxiv.org/abs/2505.02394v1) | ≈6.4908e-7<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈5.9377e-7<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim47/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 8.52% | 8.83% |
| 48<br>**L:** [≈2.3178e-8](https://doi.org/10.4153/CJM-1971-081-3) | ≈4.5291e-7<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈4.0691e-7<br>[Reoptimized nonnegative Gegenbauer cyclic kernel + signed logarithmic profile + retained Fourier witness](verifications/general_angular_20261117/) · [Sources & formulas](COMPARISONS.md) | 10.15% | 10.70% |

The 4D comparison uses the pinned Felpix-Studios certificate and its complete full-mode replay; extra digits show the small strict improvement. The 16D strict bound is 0.024900975, compared with the paper’s displayed 0.02492121. The full widened-shell count certificate is included and independently replayed. No worldwide-priority claim is made. Other densities are shown approximately to five significant digits and percentages to two decimals. Rounded display values are not substitute rigorous bounds. Full values and certificate bounds remain in [RESULTS.json](RESULTS.json); percentages use those unshortened values. The 5D, 6D, 7D and 9–16D baselines are pinned paper-reported decimal values, not independently replayed here.

These certificates use the existing [Cohn–de Laat–Salmon continuum framework](https://arxiv.org/abs/2206.15373). Proof notes state the analytic hypotheses and normalization; exact rational witnesses and independent interval verifiers check the computational conditions over full domains. Historical numerical comparisons in package notes are pinned comparisons, not current best-known claims.

<!-- lattice-lower-1-96:start -->
## Lattice-only lower bounds in 1D–96D

Verified improvements over the named prior lattice comparisons. Densities are
occupied-volume fractions. Periodic packings that are not lattices are excluded.
[Full 1D–48D catalog](LATTICE_LOWER_BOUNDS_1_48.md) ·
[Full 49D–96D catalog](LATTICE_LOWER_BOUNDS_49_96.md)

| n | Prior lattice density / source | Our lattice density / construction | Increase |
|---:|---|---|---:|
| 29 | ≈2.78789874568e-5<br>[B29 lattice — Bacher (1996)](https://doi.org/10.1007/S002220050180) · [Cohn’s lattice table](https://cohn.mit.edu/sphere-packing/#lattice-data) | **≥3.4144646907422e-5**<br>[Quaternary-code lattice section, exact center density `1/sqrt(2)`](research/quaternary_29_30_20261122/) | **22.47%** |
| 30 | ≈1.44286438425e-5<br>[Q30 lattice — Quebbemann (1984)](https://doi.org/10.1112/S0025579300010731) · [Cohn’s lattice table](https://cohn.mit.edu/sphere-packing/#lattice-data) | **≥2.1915353447830e-5**<br>[Quaternary-code lattice, exact center density `1`](research/quaternary_29_30_20261122/) | **51.89%** |
| 34 | ≈1.76697388915406e-6<br>[Cohn’s pinned table, ref. 13](https://cohn.mit.edu/sphere-packing/#13) | **≥1.83644820734493e-6**<br>[Hermitian E7 glue: primitive section; center density `4/sqrt(3)`](research/e7_hermitian35_20261123/) | **3.93%** |
| 35 | ≈9.46190415115372e-7<br>[Cohn’s pinned table, ref. 13](https://cohn.mit.edu/sphere-packing/#13) | **≥1.33811531764358e-6**<br>[Five-block Hermitian E7 glue; center density `4`](research/e7_hermitian35_20261123/) | **41.42%** |
| 36 | ≈7.33034778316707e-7<br>[Previous 36D metric-deformation lattice](constructions/dim36/) (listed by Cohn) | **≥1.11031619697710e-6**<br>[Hermitian E7 glue: order-four lift; center density `8`](research/e7_hermitian35_20261123/) | **51.47%** |
| 37 | ≈3.21356200759300e-7<br>[Cohn’s pinned table, ref. 13](https://cohn.mit.edu/sphere-packing/#13) | **≥5.24772478359786e-7**<br>[Primitive 37D section of the coupled E7 lift; center density `16/sqrt(3)`](research/e7_37_section_20261123/) | **63.30%** |
| 38 | ≈2.17566730891289e-7<br>[Cohn’s pinned lattice table, ref. 15](https://cohn.mit.edu/sphere-packing/#15) | **≥2.59631836184760e-7**<br>[Coupled E7 lift; center density `8*sqrt(2)`](research/e7_lattice38_20261123/) | **19.33%** |
| 49 | ≈4.12883717059e-9<br>[Thompson’s P48q x classical 1D lattice](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#dates-for-derived-comparators) (derived; 1971/1982 sources) | **≥7.15135575565e-9**<br>[P48q half-coset lift](research/covering_and_deep_cosets_20261115/DEEPER_PROOF.md) | **+73.20%** |
| 50 | ≈8.40815545965e-10<br>[Thompson’s P48q x classical 2D lattice](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#dates-for-derived-comparators) (derived; 1971/1982 sources) | **≥2.52244663789e-9**<br>[P48q deep plane with fitted vertical metric](research/adaptive_layers_20261115/results.json) | **+200.00%** |
| 51 | ≈1.79842551173e-10<br>[Thompson’s P48q x classical 3D lattice](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#dates-for-derived-comparators) (derived; 1971/1982 sources) | **≥9.34489307985e-10**<br>[P48n quaternionic extension](research/p48n_hurwitz_20261122/README.md) | **+419.61%** |
| 52 | ≈1.09497963768e-10<br>[Leterrier54, Corollary 1.4 / Blichfeldt section to 52D](https://doi.org/10.1007/s40993-022-00321-0) (derived; 2022 source) | **≥3.95931817049e-10**<br>[P48n quaternionic extension](research/p48n_hurwitz_20261122/README.md) | **+261.58%** |
| 53 | ≈1.25612580964e-10<br>[Leterrier54, Corollary 1.4 / Blichfeldt section to 53D](https://doi.org/10.1007/s40993-022-00321-0) (derived; 2022 source) | **≥1.66176494990e-10**<br>[P48n quaternionic extension](research/p48n_hurwitz_20261122/README.md) | **+32.29%** |
| 57 | ≈8.85119018167e-12<br>[Shioda MW56 x classical 1D lattice](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) (derived; 1991/2019 sources) | **≥1.61599884109e-11**<br>[MW56 depth-7 half-coset lift](research/packing_binary57_58_20261117/README.md) | **+82.57%** |
| 58 | ≈4.47493552650e-12<br>[Classical Gamma72/Leech-block section](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) (derived; 2010/2012 sources) | **≥5.58243739609e-12**<br>[MW56 depth-7 complex-multiplication lift](research/packing_binary57_58_20261117/README.md) | **+24.74%** |
| 60 | [1.66896497767e-12, 1.68057354713e-12]<br>[Elkies E60, reported by Vardy (1999), §4](https://doi.org/10.1090/S0002-9947-99-02169-8) | **≥1.92394316159e-12**<br>[Γ72 section normal to √2 K12](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#60d-coxetertodd-section) (derived; 1953/2010/2012 sources) | **+14.48%** |
| 65 | ≈3.20137853771e-13<br>[Classical Gamma72/Leech-block section](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) (derived; 2010/2012 sources) | **≥3.89626482463e-13**<br>[Ne64 deep order-three lift](research/novel_routes_20261113/higher_parent/ne65_verified.json) | **+21.70%** |
| 73 | ≈2.13231497472e-14<br>[Nebe Gamma72 x classical 1D lattice](https://arxiv.org/abs/1008.2862) (derived; 2010/2012 sources) | **≥4.26462994945e-14**<br>[Γ72 deep half-coset lift](research/gamma72_depth6_20261122/README.md) | **+100.00%** |
| 74 | ≈3.57517773793e-15<br>[Nebe Gamma72 x classical 2D lattice](https://arxiv.org/abs/1008.2862) (derived; 2010/2012 sources) | **≥1.43007109517e-14**<br>[Γ72 deep half-coset lift](research/gamma72_depth6_20261122/README.md) | **+300.00%** |
| 75 | ≈6.31575054665e-16<br>[Nebe Gamma72 x classical 3D lattice](https://arxiv.org/abs/1008.2862) (derived; 2010/2012 sources) | **≥5.05260043732e-15**<br>[Γ72 deep half-coset lift](research/gamma72_depth6_20261122/README.md) | **+700.00%** |
| 76 | ≈1.27986613975e-16<br>[Nebe Gamma72 x classical 4D lattice](https://arxiv.org/abs/1008.2862) (derived; 2010/2012 sources) | **≥2.04778582360e-15**<br>[Γ72 deep half-coset lift](research/gamma72_depth6_20261122/README.md) | **+1500.00%** |
| 77 | ≈1.03072807892e-16<br>[Bachoc–Nebe L80 section normal to 2A3](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#77d-bachocnebe-section) (derived; 1998 source) | **≥6.73268761846e-16**<br>[Half-root Γ72–D4 tilt](research/gamma72_root_tilt_20261122/README.md) | **+553.19%** |
| 78 | ≈9.52377596140e-17<br>[Bachoc--Nebe L80 bare section, normal rank 2](research/cross_dimension_lower_20261115/l80_results.json) (derived; 1998 source) | **≥2.53967358970e-16**<br>[Half-root Γ72–D4 tilt](research/gamma72_root_tilt_20261122/README.md) | **+166.66%** |
| 79 | ≈9.27474267490e-17<br>[Bachoc--Nebe L80 perpendicular to a minimum vector](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf) (derived; 1998 source) | **≥1.10606606099e-16**<br>[MW80 primitive section normal to 162/41](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#79d-bare-mw80-section) | **+19.25%** |
| 81 | ≈1.59531159360e-17<br>[Shioda MW80 x classical 1D lattice](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) (derived; 1991/2019 sources) | **≥2.98455470421e-17**<br>[MW80 half-coset lift](research/mw80_half_coset_81_20261123/README.md) | **+87.08%** |
| 82 | ≈2.54181457849e-18<br>[Shioda MW80 x classical 2D lattice](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) (derived; 1991/2019 sources) | **≥9.17033181889e-18**<br>[A2 section of the 83D extension](research/l80_a2_section_20261122/README.md) | **+260.77%** |
| 83 | ≈5.44626155723e-19<br>[BCS88 section normal to √2 D5](research/bcs83_comparator_20261117/README.md) (derived; 1982 source) | **≥3.08087078338e-18**<br>[Bachoc–Nebe parity extension](research/packing83_85_deep_verified_20261117/README_DEEP.md) | **+465.68%** |
| 84 | ≈2.97020393266e-19<br>[BCS88 root-normal section (4 normals)](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) (derived; 1982 source) | **≥1.18808157306e-18**<br>[Bachoc–Nebe parity extension](research/packing83_85_deep_verified_20261117/README_DEEP.md) | **+300.00%** |
| 85 | ≈1.61034685155e-19<br>[BCS88 root-normal section (3 normals)](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) (derived; 1982 source) | **≥4.55474871517e-19**<br>[Bachoc–Nebe parity extension](research/packing83_85_deep_verified_20261117/README_DEEP.md) | **+182.84%** |
| 86 | ≈1.00229841948e-19<br>[BCS88 root-normal section (2 normals)](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) (derived; 1982 source) | **≥2.61014583751e-19**<br>[Primitive MW92 kernel section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | **+160.41%** |
| 87 | ≈6.57892648392e-20<br>[BCS88 root-normal section (1 normals)](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) (derived; 1982 source) | **≥2.04757003794e-19**<br>[Primitive MW92 kernel section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | **+211.23%** |
| 88 | ≈4.95809518917e-20<br>[Bos–Conway–Sloane BCS88 (1982)](https://neilsloane.com/doc/Me103.pdf) | **≥1.85638056287e-19**<br>[Primitive MW92 kernel section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | **+274.41%** |
| 89 | [5.33520442881e-20, 5.33661675466e-20]<br>[Guaranteed MW92 section via Rogers’ inequality](research/covering_and_deep_cosets_20261115/SECTIONS_PROOF.md) (derived; 1958/1991/2019 sources) | **≥1.52239236465e-19**<br>[Primitive MW92 kernel section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | **+185.27%** |
| 90 | [7.10853349171e-20, 7.10978166244e-20]<br>[Guaranteed MW92 section via Rogers’ inequality](research/covering_and_deep_cosets_20261115/SECTIONS_PROOF.md) (derived; 1958/1991/2019 sources) | **≥1.43364905437e-19**<br>[Primitive MW92 kernel section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | **+101.64%** |
| 91 | [9.59862180378e-20, 9.59946031324e-20]<br>[Guaranteed MW92 section via Rogers’ inequality](research/covering_and_deep_cosets_20261115/SECTIONS_PROOF.md) (derived; 1958/1991/2019 sources) | **≥1.28046110298e-19**<br>[Primitive MW92 kernel section](research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels) | **+33.38%** |
| 93 | ≈1.70235478236e-20<br>[Shioda MW92 x classical 1D lattice](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) (derived; 1991/2019 sources) | **≥3.045564894259504e-20**<br>[MW92 degree-13 covering refinement; lattice existence proof](research/mw92_cfr_taylor13_20261123/) | **+78.9030%** |
| 94 | ≈2.53431409158e-21<br>[Shioda MW92 x classical 2D lattice](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) (derived; 1991/2019 sources) | **≥8.036027386567117e-21**<br>[MW92 degree-13 covering refinement; lattice existence proof](research/mw92_cfr_taylor13_20261123/) | **+217.0888%** |
| 95 | ≈3.98071737750e-22<br>[Shioda MW92 x classical 3D lattice](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) (derived; 1991/2019 sources) | **≥2.104056220143222e-21**<br>[MW92 degree-13 covering refinement; lattice existence proof](research/mw92_cfr_taylor13_20261123/) | **+428.5620%** |
| 96 | ≈2.79502329482e-22<br>[BCS96, Table 1 (1982)](https://neilsloane.com/doc/Me103.pdf) | **≥5.93710931559e-22**<br>[MW92 cusp reduction glued to D4](research/mw92_cusp_d4_20261122/README.md) | **+112.41%** |

The 29D and 30D rows improve Cohn’s listed lattice bounds and tie his listed
unrestricted bounds. The 34–37D E7 rows improve both pinned tables. The 36D
gain is measured against the preceding Cohn-listed metric-deformation lattice.
The 38D lattice matches the previously verified periodic density. Its 37D section improves the prior 37D density.
From 49D onward, the gains are rounded down and use the upper end of a prior
comparison interval where applicable. Years in derived rows date the cited
sources, not publication of the derived bound. The catalogs include ties, source
notes and exact center densities. The comparisons do not establish global
optimality or exhaustive historical priority.
<!-- lattice-lower-1-96:end -->

<!-- lattice-upper-20261116:start -->
## Lattice-only upper bounds

The current register gives certified upper bounds in **9–48D** and compares them with Zubrilina’s **13–48D** repository and our prior registry. The latest high-dimensional bounds jointly optimize the Fourier test and the interactions of short vectors in each lattice plane. In 13D a determinant-denominator bound gives a strict improvement of about 8.03 × 10⁻⁶³⁶% over the exact Mordell value. In 25D, 31D and 32D the Fourier test also bounds how many short vectors can occupy selected squared-length intervals. [Short-vector count certificates](verifications/lattice/poisson_anchor_windows_20261122/) · [Joint certificate and replay](verifications/lattice/joint_fourier_plane_20261122/) · [13D proof](verifications/lattice/perfect_determinant_13/) · [Full comparison](LATTICE_RESULTS.md#current-lattice-upper-bounds-in-dimensions-948) · [Machine-readable bounds](LATTICE_UPPER_BOUNDS_9_48.json)

| Dimension | Current lattice density upper | Below Zubrilina’s 13–48D bound |
|---:|---:|---:|
| 14 | **0.033996852420018639356** | 0.00008% |
| 15 | **0.025163909222702077102** | 0.000026% |
| 16 | **0.017231432341539144664** | 2.89e−10% |
| 19 | **0.0076312827957442716237** | 8.3193e−12% |
| 23 | **0.0025265652019047523943** | 1.9131e−12% |
| 25 | **0.0012865460837431921899** | 0.00093% |
| 31 | **0.00013541049750198333762** | 0.01411% |
| 32 | **0.000095866360403730692796** | 0.010993% |
| 36 | **0.000027145852648811524495** | 0.000109% |
| 37 | **0.00001817642454038903306** | 0.378933% |
| 38 | **0.00001229554573503757040** | 0.383695% |
| 39 | **0.00000845291367154137428** | 0.282223% |
| 42 | **0.00000346603735498588186** | 0.18193% |
| 48 | **0.00000034510710122965769** | 0.057497% |

The exact comparison gives **35 smaller upper bounds** and **one tie**, the exact 24D optimum, across the 36 benchmark dimensions. The largest reduction is **0.383695% in 38D**. [Pinned benchmark](https://github.com/AstroNinja389/sphere-packings/blob/7b63ee2094b7489aea3da83349c547c8da4cb388/proofs/COMPARISON.md) · [Historical project registry](LATTICE_RESULTS.md#historical-project-registry)
<!-- lattice-upper-20261116:end -->

## Recent updates

**22 November:** a [quaternary-code lattice and its section](research/quaternary_29_30_20261122/) match the unrestricted lower-bound densities in **29D** and **30D**, improving the immediately preceding lattice bounds by **1.55%** and **28.09%**, respectively.

**22 November:** shifted Bacher lattices improve the catalogued **29D and 30D lattice-only density lower bounds by 20.60% and 18.58%**, respectively. The unrestricted B29* and T30 bounds are unchanged. [Exact construction and replay](research/bacher_slip_29_30_20261122/).
**22 November:** the summary table now reflects the already verified **82D density lower bound of 9.17033181889404e-18**, 7.0436% above the stale mixed-depth entry. The [principal-section proof](research/l80_a2_section_20261122/) is a simpler alternative to the [earlier all-deep binary-plane proof](research/packing83_85_deep_verified_20261117/README.md), not a new density record.

**22 November:** explicit lattices based on Nebe’s Γ72 improve the previous lower bounds in **73–77D** by **9.5445%, 20%, 31.4534%, 44% and 33.3180%**, respectively. The same construction gives a new 78D lattice lower bound of **1.39319236921112e-16**. [Construction and exact replay](research/gamma72_leech_layers_20261122/).

**22 November:** explicit minimum-6 lattices improve the preceding **51–53D** lower bounds by **17.2603%**, **37.5%** and **11.0281%**, respectively. A nine-translate periodic packing improves **96D** by **2.2869%**. New lattice existence bounds cover **93–96D**. [51–53D proof](research/p48n_hurwitz_20261122/) · [93–96D proofs](research/mw92_extensions_20261122/) · [Complete 1–96D lower-bound table](LOWER_BOUNDS_1_96.md).

**18 November:** seven independently verified unrestricted upper-bound improvements in **17–23D**, reducing the previous repository upper bounds by **0.319–1.498%**. [Exact before/after values and verification scope](CHANGES_20261118.md).

<details>
<summary><strong>Week of 16 November 2026</strong> · stronger upper bounds and new packing constructions</summary>

**17 November:** unrestricted 49–54D certificates were followed by stronger compact-kernel refinements, while three-feature kernels and dimensional transfers tightened lattice-only bounds. Packing work added advances in 57–58D, 81–90D and 94–95D; explicit constructions and existence results remain distinguished in their tables. [Chronological release notes and exact comparisons](CHANGES_20261117.md).

Latest verified values from this week:

| Scope | Current release | Proof and comparison |
| --- | --- | --- |
| Unrestricted upper, 58D | **<1.073633e-8** | [Complete audit](verifications/general58_pc29_geomfix_20261117/); [numerical comparison](NUMERICAL_BENCHMARKS49_96.md#57d-and-58d-certified-bounds-with-numerical-literature-comparisons), no SOTA claim |
| Unrestricted upper, 57D | **<1.5718326e-8** | [Complete audit](verifications/general57_pc29_geomfix_20261117/); [numerical comparison](NUMERICAL_BENCHMARKS49_96.md#57d-and-58d-certified-bounds-with-numerical-literature-comparisons), no SOTA claim |
| Unrestricted upper, 54D | **<4.8695034e-8** | [Complete audit](verifications/general54_pc29_geomfix_20261117/); [numerical benchmark comparison](NUMERICAL_BENCHMARKS49_96.md) |
| Unrestricted upper, 53D | **<6.8841195e-8** | [Complete audit](verifications/general53_pc29_geomfix_20261117/); 7.0662% below the previous repository certificate |
| Unrestricted upper, 51D | **<1.37012082e-7** | [Joint scalar/kernel and complete audit](verifications/general51_jointp6_20261117/); 0.5584% below the preceding certificate |
| Unrestricted upper, 52D | **<9.6610185e-8** | [Three-feature kernel and complete audit](verifications/general52_rank3mid_20261117/); 0.8519% below the preceding certificate |
| Unrestricted upper, 50D | **<1.99548924e-7** | [Complete audit](verifications/general50_pc29_20261117/) |
| Lattice-only upper, 36D | **<2.727705e-5** | [Dimension restriction and free-origin optimization](verifications/lattice/compact_restriction36_20261117/); 0.457675% below the preceding certificate |
| Lattice-only upper, 40–42D | 40D: **<5.9777359372e-6**<br>41D: **<4.3927514316e-6**<br>42D: **<3.473037643e-6** | [Dimension-specific compact kernels](verifications/lattice/compact_power6_40_42_20261117/); complete independent audits |
| Lattice-only upper, 37–39D | 37D: **<1.8268708324e-5**<br>38D: **<1.23521549304e-5**<br>39D: **<8.4799401622e-6** | [Dimension-specific compact kernels](verifications/lattice/compact_power6_37_39_20261117/); complete independent audits |
| Lattice-only upper, 43–44D | 43D: **<2.5363822051e-6**<br>44D: **<1.7676809433e-6** | [Dimension-specific three-feature kernels](verifications/lattice/compact_power6_adapted_20261117/); complete independent audits |
| Lattice-only upper, 45–48D | 45D: **<1.1728551765e-6**<br>46D: **<7.721596228e-7**<br>47D: **<5.130033534e-7**<br>48D: **<3.4559190288e-7** | [Three-feature kernels and exact dimensional transfers](verifications/lattice/compact_power6_45_48_20261117/) |
| Packing lower, 94D/95D | **≥7.84434540405174e-21; ≥2.06770138493108e-21** | [Anisotropic torus-offset flags](research/mw92_flags9495_verified_20261117/); periodic existence, not explicit offsets or lattice claims |

[49–96D numerical benchmark comparisons](NUMERICAL_BENCHMARKS49_96.md) are explicitly **not certified prior upper bounds**. [Release history](RELEASE_NOTES_20261117.md) preserves earlier announcements, superseded values and proof links; [17 November changes](CHANGES_20261117.md) records exact updates.

</details>

<details>
<summary><strong>Week of 9 November 2026</strong> · initial constructions and successive upper-bound refinements</summary>

- **12 November:** added the explicit 36D packing improvement over KP36 and upper-bound certificates in 9D, 10D and 16D. [Original releases and proofs](RELEASE_NOTES_20261117.md).
- **14 November:** strengthened upper bounds in 16 dimensions, using signed logarithmic profiles and radial enrichment. [Exact before/after values](CHANGES_20261114.md).
- **15 November:** extended the attributed 4D certificate, added further upper-bound refinements, and verified an explicit 59D packing with a stronger literature-derived comparison. [Chronology, attribution and replay](CHANGES_20261115.md).

Proof dependencies remain in each package. In particular, the 9D certificate retains external average-neighbor bound 373; the 36D minimum proof uses outward-rounded binary64 enumeration under documented compiler/arithmetic assumptions, with an independent rerun of the same audited implementation. The 4D result extends the attributed pinned Felpix-Studios certificate. Packing constructions are lower bounds, not optimality claims.

</details>

[Full release history and archived research notes](RELEASE_NOTES_20261117.md). Historical values may be superseded; the current result tables and linked certificates govern.

## Reproduce

For the upper-bound suite, Python 3.12 or newer and a C++20 compiler are required by the pinned 4D proof. Install dependencies and run:

    python -m pip install -r requirements.txt
    python run_all.py

Never disable assertions (`python -O`). No optimization solver or network is needed after dependency installation. The 48D compact triangle replay may take several minutes. Logs record the original audits; fresh replay is recommended. No GitHub Actions workflow is installed.

For the latest 34–37D lower bounds, run `python research/e7_hermitian35_20261123/verify.py` and `python research/e7_hermitian35_20261123/basis.py`. The separate historical 36D lower-bound replay is `python constructions/dim36/verify_all.py`; it requires SymPy and a C++17 compiler and performs the complete minimum enumeration. It is not silently included in the historical upper-bound suite.

The [51–53D lattice certificate](research/p48n_hurwitz_20261122/#replay) includes a full replay of its coset exclusions. The [93–96D proofs](research/mw92_extensions_20261122/#replay) and the [1–96D table](LOWER_BOUNDS_1_96.md) can be checked with the Python standard library.

## Layout

- [verifications/](verifications/): unrestricted and lattice-only upper-bound certificates, with package-specific proofs and replay instructions.
- [RESULTS.json](RESULTS.json) and [LATTICE_RESULTS.json](LATTICE_RESULTS.json): separate machine-readable registries. [Additional high-dimensional releases](NUMERICAL_BENCHMARKS49_96.md#certified-results-compared-with-these-numerical-benchmarks) link their self-contained certificates.
- [constructions/](constructions/) and [research/](research/): packing constructions, source audits and experiments. A research directory alone is not a verified result; consult its stated status and proof.
- [Root replay suite](run_all.py) and [lattice replay suite](verifications/lattice/replay.py): separate entry points. Later releases may have their own runners; follow each package’s README rather than assuming the root suite covers everything.
- [Release history](RELEASE_NOTES_20261117.md): chronology, superseded witnesses and detailed audit notes.

The referenced continuum theorem and standard positive-definite Gegenbauer addition theorem are mathematical inputs, not independently formalized theorems in this repository.

[Origin-coordinate](verifications/origin_slack/) and [cyclic-origin refinements](verifications/cyclic_origin_slack/) remain available as earlier proof stages. Use the registries and current release links for the latest certificates.

## Research tools

The [research index](research/README.md) separates registered certificates,
historical results and unsuccessful searches. The [16 November branch
reconciliation](research/branch_integration_20261116/README.md) records where
the three research branches were incorporated and which results were superseded.

- [Dimension-generic annular three-point search](research/annular_three_point/): global-positivity SDP discovery, exact export and full-domain checks. The degree-31 18D witness is now registered after independent audit. Degree-23 17D/18D and degree-31 17D certificates are weaker than published bounds; 22D/23D candidates are rejected. All earlier experiments remain preserved.
- [Origin-coordinate slack research](research/origin_slack/): exact Schur-complement postprocessing and independent 5D audit. This preserves the initial research snapshot; the subsequent fully replayed refinements are registered in [verifications/origin_slack](verifications/origin_slack/).

Research dependencies are separate from the verification-only root requirements. See each package for reproducibility commands, tested versions, scope and checksums.

The newly audited [18D certificate](verifications/dim18/) proves density <0.013408931, improving the pinned published upper by 0.566611%. It uses the continuum theorem without external neighbor-count bounds. Comparison freshness and worldwide-priority caveats remain in force.

The [17D solid-harmonic cyclic certificate](verifications/dim17/) proves density <0.018335552. Its proof uses an explicit cyclic generalization and periodic PSD summation, not the individual-sign hypothesis of the printed Theorem 1.1. The [search code and rejected experiments](research/solid_cyclic/) are retained separately.

The [18D cyclic refinement](verifications/dim18/refinement_cyclic/) strengthens the earlier bound to <0.013284728. Both cyclic packages retain explicit summation proofs and independent exact-domain audits; earlier witnesses remain available.

The [19D cyclic refinement](verifications/dim19/refinement_cyclic/) proves density <0.009557491 and the [20D cyclic refinement](verifications/dim20/refinement_cyclic/) proves <0.006977702. Both use the explicitly proved cyclic generalization and independently reconstructed exact-domain checks. All older witnesses and Schur certificates remain replayable.

The current 17D bound <0.018204096 uses a [mixed-Gaussian repaired cyclic witness and separate saturation](verifications/dim17/refinement_repaired/). The [common cyclic origin saturation](verifications/cyclic_origin_slack/) sharpens 18D, 19D and 20D to <0.013224473, <0.009521000 and <0.006936042. All prior witnesses remain preserved and replayed. See [methods and experiments](METHODS_AND_EXPERIMENTS.md) for the systematic distinction between valid proof transformations and unverified search directions.

The [47D restriction and origin-Gram certificate](verifications/dim47/) gives <6.16879891e-7 by transferring the fully replayed 48D witness. The pinned upper decreases 4.96%, and the additive gap to the cited 2025 antipode construction decreases 5.14%. The proof and independent verifier explicitly preserve the cyclic theorem and Fourier restriction hypotheses.

The [6D count-enhanced certificate](verifications/dim6/) gives <0.410189987, using externally published near-neighbor bounds. The [45D](verifications/dim45/) and [46D](verifications/dim46/) exact radial/angular certificates give <1.311423575e-6 and <9.14882847e-7. All new proof chains and preserved historical packages are included in the aggregate replay.

The [portable 6D discovery/export package](research/dim6_discovery/) reproduces the successful bounded warm-start search and exports its frozen candidate for full exact verification.

The [improved 6D certificate](verifications/dim6/refinement_mixture99/) strengthens the earlier result to <0.409135399, closing 3.09% of the pinned additive gap. Its physical and Fourier signs hold without Gaussian repair; the same published neighbor-count dependencies remain. The earlier 6D proof is preserved.

The [6D refinement discovery package](research/dim6_refinement_discovery/) reconstructs the frozen 99/100 source exactly from preserved checkpoints and exports it for full verification.

The [44D Gaussian cyclic certificate](verifications/dim44/) gives <1.842868599e-6, about 3.10% below the pinned upper. The [small 6D follow-on](verifications/dim6/refinement_patch5/) gives <0.409128533, with no Gaussian repair and unchanged external count dependencies.

The [7D count-enhanced certificate](verifications/dim7/) gives <0.320584772, closing 2.05% of the pinned additive gap to E7. Its external near-neighbor count inputs are cited explicitly.

The [guarded 5D refinement](verifications/dim5/refinement_guarded/) gives <0.510129239; its unguarded checkpoint is invalid. The [22D tetrahedral aggregate research](research/tetrahedral_aggregate/) proves stronger auxiliary edge/triangle inequalities but makes no density-record claim.

Portable exact-source reconstruction and export companions are available for [guarded 5D](research/dim5_guarded_discovery/), [7D](research/dim7_discovery/) and the [small 6D follow-on](research/dim6_patch5_discovery/). Numerical checkpoints remain distinct from certified witnesses.

The [5D long-edge lens extension](verifications/dim5/refinement_lens/) gives <0.509996147, closing 4.83% of the pinned additive gap. It uses newly audited long-edge PSD accounting, not the old pointwise far-sign condition. Earlier proof packages remain preserved.

The [stronger 5D spherical-lens follow-on](verifications/dim5/refinement_spherical_lens/) gives <0.509953498, closing 4.92% of the pinned additive gap. Its exact transverse spherical-cap proofs and new tail budgets are replayed; the centroid-only version remains available. Portable search/export companions are included for [centroid](research/dim5_lens_discovery/) and [spherical](research/dim5_spherical_lens_discovery/) lenses.

The [22D](verifications/dim22/), [23D](verifications/dim23/), [26D](verifications/dim26/) and [27D](verifications/dim27/) Gaussian cyclic certificates have fresh primary and independent proof replays. [Release checks](research/gaussian_22_27_release/) distinguish those from reused historical audits. Portable discovery and exact export are included for [22/23](research/gaussian22_23_discovery/), [26](research/dim26_gaussian_discovery/) and [27](research/dim27_gaussian_discovery/).

The remaining [29D](verifications/dim29/), [30D](verifications/dim30/), [31D](verifications/dim31/), [33D](verifications/dim33/), [35D](verifications/dim35/) and [37D](verifications/dim37/) entries complete the Gaussian cyclic certificates in dimensions 25–40. Their [release checks](research/gaussian_remaining_release/) distinguish fresh primary/export replays, complete independent audits and reused unchanged historical proofs.

The [21D Gaussian cyclic refinement](verifications/dim21/refinement_gaussian/) proves density<0.004985977000, retaining all earlier21D proofs. Its [release evidence](research/gaussian21_release/) distinguishes fresh proof checks from unchanged historical audits.

The [21D discovery companion](research/dim21_gaussian_discovery/) includes the original seven-solve search history, exact inputs and deterministic export replay; numerical checkpoints are not proof certificates.

The [46D](verifications/dim46/refinement_angular/) and [47D](verifications/dim47/refinement_angular/) angular refinements lower the origin cost with rebuilt target-dimensional angular kernels and unchanged radial witnesses. Earlier proofs are preserved; [release checks](research/angular_border_release/) document fresh new and historical46/47 replays and reuse of other unchanged audits.

The subsequent [46D](verifications/dim46/refinement_joint/) and [47D](verifications/dim47/refinement_joint/) joint refits reoptimize the radial Fourier pair for the rebuilt angular kernel, improving both bounds further. [Release evidence](research/angular_radial_refit_release/) includes fresh primary/independent replays, exact export checks and explicit historical-check scope.

## Citation

If you use these bounds, packing constructions, or verification packages, please cite:

Alexey Kravatskiy. (2026). *Sphere-packing density: lower and upper bounds for unrestricted and lattice packings*. https://github.com/alexlegeartis/SpherePackings

[CITATION.cff](CITATION.cff) supplies GitHub's citation metadata; [CITATION.bib](CITATION.bib) provides BibTeX. For a reproducible reference to particular results, also record the commit used.


## AI-assisted development

The constructions, the verification code and the write-ups in this repository were developed by the author working with Claude Code (Anthropic), ChatGPT Pro (GPT 6), and Astra & Sol 6.1 in Codex.

