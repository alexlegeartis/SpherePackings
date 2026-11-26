# The 22-dimensional Conway–Sloane projection: scope and exact limits

No new density record. We use the [22D Conway–Sloane comparator](https://arxiv.org/pdf/math/0207182) with center density `3/sqrt(12) * (11/12)^11`.

- [Anticode theorem](ANTICODE.md): in the fixed projected Leech grid `Λ22*`, every subset, periodic or not, with squared minimum at most 4 has center density at most the comparator. The restriction to this carrier and this squared minimum is essential.
- [Periodic equality](PERIODIC_EQUALITY.md): at the same maximum a full-rank-periodic subset is the original three-coset packing up to translation and inversion. The forcing set uses 23 exact four-point stencils; no assertion for aperiodic equality.
- [General torus diameter guarantee](TORUS_DIAMETER.md): Conway–Sloane attains the optimal guarantee from ambient minimum and normal diameter over all rank-two Leech sections and arbitrary numbers of labels, including winding. This does not bound the actual horizontal metrics when a shell is missing.
- [Local rigidity in the original three-translate cell](LOCAL_RIGIDITY.md): strict local density maximum under joint metric deformation and independent layer shifts, without larger supercell claims.
- The separate search over 430 tested doubled finite phase palettes is exactly negative. This is a test of those 430 palettes, not every carrier or all off-grid shifts. `rows.json.gz` is the exact candidate input, `witnesses.json` supplies 603 explicit fellow-vector corrections for the 90 cases not closed by the first shell. The replay verifies all originally tested palettes and every decisive clique threshold.

```bash
python3 verify_anticode.py
python3 verify_periodic_equality.py
node verify_torus.js
node verify_local_rigidity.js
node verify_double_palettes.js rows.json.gz witnesses.json
```

The first four programs are self-contained (Python standard library or Node). The final program regenerates the Leech minimal shell, reads compressed rows with Node's standard zlib, and uses only integer arithmetic for all decisive comparisons. An optional discovery run is in `scan_double_palettes.py` (NumPy); it is not needed for verification.
