# Lattice sphere-packing lower bounds, dimensions 1–48

[Companion lattice-only catalog for 49D–96D](LATTICE_LOWER_BOUNDS_49_96.md).

These are densities attained by **lattice packings**. A genuinely nonlattice
union of translates does not qualify. The baseline is the lattice-only table on
[Henry Cohn’s sphere-packing page](https://cohn.mit.edu/sphere-packing/#lattice-data),
checked on 22 November 2026. The [unrestricted lower-bound catalog](LOWER_BOUNDS_1_96.md)
is separate. This table does not claim optimality or exhaustive priority.

## Updates for Cohn’s lattice table

| Dimension | Listed center density | New exact center density | Occupied density (lower bound) | Gain vs. listed lattice | Construction |
|---:|---:|---:|---:|---:|---|
| 29 | `3^(-1/2)` | **`1/sqrt(2)`** | **3.4144646907422e-5** | **22.47%** | [Section of quaternary-code lattice](research/quaternary_29_30_20261122/) |
| 30 | `2^-22 * 3^(27/2)` | **`1`** | **2.1915353447830e-5** | **51.89%** | [Quaternary-code lattice](research/quaternary_29_30_20261122/) |

The 29D and 30D lattices are strictly denser than B29 and Q30, respectively,
and match the listed unrestricted B29* and T30 lower bounds. They supersede
the repository’s [shifted Bacher lattices](research/bacher_slip_29_30_20261122/)
by 1.55% and 28.09%, respectively.

The repository’s [36D lattice](constructions/dim36/) improved on KP36 by
18.97%; Cohn already lists that result in both tables, so it is not a new
update for his table.

## Full 1–48D catalog

For a lattice with shortest nonzero vector of length `2r` and Gram matrix
`G`, the center density is `r^n/sqrt(det(G))`. The occupied density is the
center density times the unit-ball volume `pi^(n/2)/Gamma(1+n/2)`. Decimal
occupied densities below are rounded downward. Source numbers point to the
construction references in the pinned Cohn table, as listed on the review date.

| Dimension | Exact center density | Occupied density (lower bound) | Construction or source |
|---:|---:|---:|---|
| 1 | `2^-1` | 1 | Integer lattice |
| 2 | `2^-1 * 3^(-1/2)` | 9.0689968211710e-1 | [Cohn, ref. 26](https://cohn.mit.edu/sphere-packing/#26) |
| 3 | `2^(-5/2)` | 7.4048048969306e-1 | [Cohn, ref. 18](https://cohn.mit.edu/sphere-packing/#18) |
| 4 | `2^-3` | 6.1685027506808e-1 | [Cohn, ref. 22](https://cohn.mit.edu/sphere-packing/#22) |
| 5 | `2^(-7/2)` | 4.6525761330925e-1 | [Cohn, ref. 23](https://cohn.mit.edu/sphere-packing/#23) |
| 6 | `2^-3 * 3^(-1/2)` | 3.7294754558206e-1 | [Cohn, ref. 23](https://cohn.mit.edu/sphere-packing/#23) |
| 7 | `2^-4` | 2.9529787314571e-1 | [Cohn, ref. 23](https://cohn.mit.edu/sphere-packing/#23) |
| 8 | `2^-4` | 2.5366950790104e-1 | [Cohn, ref. 23](https://cohn.mit.edu/sphere-packing/#23) |
| 9 | `2^(-9/2)` | 1.4577487580817e-1 | [Cohn, ref. 7](https://cohn.mit.edu/sphere-packing/#7) |
| 10 | `2^-4 * 3^(-1/2)` | 9.2021118431305e-2 | [Cohn, ref. 7](https://cohn.mit.edu/sphere-packing/#7) |
| 11 | `2^-1 * 3^(-5/2)` | 6.0432660108165e-2 | [Cohn, ref. 2](https://cohn.mit.edu/sphere-packing/#2) |
| 12 | `3^-3` | 4.9454176624244e-2 | [Cohn, ref. 14](https://cohn.mit.edu/sphere-packing/#14) |
| 13 | `2^-1 * 3^(-5/2)` | 2.9208430928107e-2 | [Cohn, ref. 29](https://cohn.mit.edu/sphere-packing/#29) |
| 14 | `2^-4 * 3^(-1/2)` | 2.1624096082447e-2 | [Cohn, ref. 27](https://cohn.mit.edu/sphere-packing/#27) |
| 15 | `2^(-9/2)` | 1.6857570656762e-2 | [Cohn, ref. 3](https://cohn.mit.edu/sphere-packing/#3) |
| 16 | `2^-4` | 1.4708164397430e-2 | [Cohn, ref. 3](https://cohn.mit.edu/sphere-packing/#3) |
| 17 | `2^-4` | 8.8113191823211e-3 | [Cohn, ref. 27](https://cohn.mit.edu/sphere-packing/#27) |
| 18 | `2^-3 * 3^(-1/2)` | 5.9283687184694e-3 | [Cohn, ref. 27](https://cohn.mit.edu/sphere-packing/#27) |
| 19 | `2^(-7/2)` | 4.1208062797686e-3 | [Cohn, ref. 27](https://cohn.mit.edu/sphere-packing/#27) |
| 20 | `2^-3` | 3.2258614237517e-3 | [Cohn, ref. 27](https://cohn.mit.edu/sphere-packing/#27) |
| 21 | `2^(-5/2)` | 2.4658847115024e-3 | [Cohn, ref. 27](https://cohn.mit.edu/sphere-packing/#27) |
| 22 | `2^-1 * 3^(-1/2)` | 2.1276601452758e-3 | [Cohn, ref. 28](https://cohn.mit.edu/sphere-packing/#28) |
| 23 | `2^-1` | 1.9053281934260e-3 | [Cohn, ref. 28](https://cohn.mit.edu/sphere-packing/#28) |
| 24 | `1` | 1.9295743094039e-3 | [Cohn, ref. 28](https://cohn.mit.edu/sphere-packing/#28) |
| 25 | `2^(-1/2)` | 6.7721200977318e-4 | [Cohn, ref. 11](https://cohn.mit.edu/sphere-packing/#11) |
| 26 | `3^(-1/2)` | 2.6922005043380e-4 | [Cohn, ref. 11](https://cohn.mit.edu/sphere-packing/#11) |
| 27 | `3^(-1/2)` | 1.2867528120269e-4 | [Cohn, ref. 1](https://cohn.mit.edu/sphere-packing/#1) |
| 28 | `2 * 3^-1` | 6.9758736616563e-5 | [Cohn, ref. 1](https://cohn.mit.edu/sphere-packing/#1) |
| 29 | **`1/sqrt(2)`** | 3.4144646907422e-5 | [Section of quaternary-code lattice](research/quaternary_29_30_20261122/) (improves Cohn’s listed lattice) |
| 30 | **`1`** | 2.1915353447830e-5 | [Quaternary-code lattice](research/quaternary_29_30_20261122/) (improves Cohn’s listed lattice) |
| 31 | `2^(-47/2) * 3^15` | 1.1837765185933e-5 | [Cohn, ref. 31](https://cohn.mit.edu/sphere-packing/#31) |
| 32 | `2^-24 * 3^16` | 1.1040749308859e-5 | [Cohn, ref. 31](https://cohn.mit.edu/sphere-packing/#31) |
| 33 | `2^-25 * 3^(33/2)` | 4.1406882896496e-6 | [Cohn, ref. 13](https://cohn.mit.edu/sphere-packing/#13) |
| 34 | `2^-25 * 3^(33/2)` | 1.7669738891540e-6 | [Cohn, ref. 13](https://cohn.mit.edu/sphere-packing/#13) |
| 35 | `2^(3/2)` | 9.4619041511537e-7 | [Cohn, ref. 13](https://cohn.mit.edu/sphere-packing/#13) |
| 36 | `16777216/3176523` | 7.3303477831670e-7 | [Metric deformation](constructions/dim36/) (already listed by Cohn) |
| 37 | `2^(5/2)` | 3.2135620075930e-7 | [Cohn, ref. 13](https://cohn.mit.edu/sphere-packing/#13) |
| 38 | `2^-19 * 3^(31/2) * 5^-1` | 2.1756673089128e-7 | [Cohn, ref. 15](https://cohn.mit.edu/sphere-packing/#15) |
| 39 | `2^-20 * 3^16 * 5^(-1/2)` | 1.6802817736543e-7 | [Cohn, ref. 15](https://cohn.mit.edu/sphere-packing/#15) |
| 40 | `2^-20 * 3^16` | 1.4798339931342e-7 | [Cohn, ref. 15](https://cohn.mit.edu/sphere-packing/#15) |
| 41 | `2^-21 * 3^17` | 8.6368303579508e-8 | [Cohn, ref. 15](https://cohn.mit.edu/sphere-packing/#15) |
| 42 | `2^-21 * 3^(35/2)` | 5.7516899055409e-8 | [Cohn, ref. 15](https://cohn.mit.edu/sphere-packing/#15) |
| 43 | `2^(-45/2) * 3^19` | 4.0157190248136e-8 | [Cohn, ref. 11](https://cohn.mit.edu/sphere-packing/#11) |
| 44 | `2^-23 * 3^20` | 3.2008535265505e-8 | [Cohn, ref. 11](https://cohn.mit.edu/sphere-packing/#11) |
| 45 | `2^(-47/2) * 3^21` | 2.5231506774470e-8 | [Cohn, ref. 11](https://cohn.mit.edu/sphere-packing/#11) |
| 46 | `2^-23 * 3^(43/2)` | 2.2717980357903e-8 | [Cohn, ref. 11](https://cohn.mit.edu/sphere-packing/#11) |
| 47 | `2^-24 * 3^23` | 2.1466078267465e-8 | [Cohn, ref. 11](https://cohn.mit.edu/sphere-packing/#11) |
| 48 | `2^-24 * 3^24` | 2.3178295310541e-8 | [Cohn, ref. 29](https://cohn.mit.edu/sphere-packing/#29) |

## Reproduce

[Source factors and reference numbers](research/lattice_lower_bounds_1_48_20261122/cohn_inputs.json)
record the lattice-only table as read on 22 November. Run
`python3 research/lattice_lower_bounds_1_48_20261122/verify_table.py` to check
this catalog (add `--write` to rebuild it). This checks density arithmetic
with rational enclosures of π; it does not re-prove the source packings.
For the new lattices, run `python3 research/quaternary_29_30_20261122/verify.py` for
the separate complete minimum and determinant proof.
