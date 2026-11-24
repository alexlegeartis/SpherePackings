# 96D periodic packing from an MW92 half phase and twelve colors

The occupied density is at least **6.367774598225213e-22**, an improvement of **1.8708%** over the [preceding repository lower bound](../mw92_cfr_taylor13_20261123/README.md). The value of this construction is exactly

    π^48 · 2^100 / (48! · 47^7 · 231).

The result is periodic, with an existential horizontal offset from a certified covering argument. It is not an explicit lattice and does not raise the lattice-only lower bound. For comparison, the [Bos–Conway–Sloane 96D lattice](https://neilsloane.com/doc/Me103.pdf) has center density 3^48/2^24. This comparison and the prior repository floor do not assert exhaustive literature priority.

Opening live SpherePackings main 294712432b1044f2e4762823297bbceb3462b498. Parent minimum squared norm 16 and covolume 47^7 are inherited from research/mw92_extensions_20261122/README.md; exact J, order-23 T, G, and rectangle E are rechecked in the archived embed.py. The covering bound and comparator are inherited from research/mw92_cfr_taylor13_20261123/README.md and adaptive_layers_20261115/COVERING_PROOF.md; tail.py independently replays them for four balls.

**Half phase.** Let h=(E[0]+E[3])/2, fixed under T with Q(2h)=322. If a horizontal nearest offset v=x-h is fixed, 2v is a nonzero vector of the fixed lattice with Gram 23 A. Gershgorin gives A≥2I, so ||v||²≥46/4>21/2. If v is not fixed, each of its 22 order-23 differences is a nonzero parent vector and has norm²≥16. Averaging gives ||v||²≥176/23, since sum ||v−T^j v||²=46 ||(I−P)v||² for the orthogonal fixed projector. Also Q(2v)=2 mod4. Therefore the only obstruction to depth²≥21/2 is y=2v of parent norm²34 or38.

**Exact rectangle.** The integral map V verified by embed.py has VV^t=188G and each block sums to zero. Each block is in RS19 modulo 47; together the two blocks modulo 4 lie in the free rank46 cyclic C4 code. These inclusions are equality by exact index: Euclidean covolume of V(L) is 188^46*47^7, of A46⊕A46 is 47, hence index 2^92*47^52, precisely the product of the RS and C4 indices. Therefore y lies in the affine class iff it has the two RS conditions and y mod8 lies in y0+2C4, with y0 patterns (1,4,5;4,6,6) on the zero, QR, NQR positions. The second block yB is even; ||yB/2||² is even (zero sum) and divisible by 47 (RS square degree≤38); hence both block energies in units 188 are even. Norm34 needs one block q≤16 (actual squared energy≤3008), norm38 one block q≤18 (≤3384). The smaller block is unique. Fixing it determines the complementary mod8 fibre; its homogeneous lattice is 4 times the appropriate full projection QR∩RS lattice.
 
**Exact enumeration (completed).** Four recorded original outer shards per side partition the first six integral GSO coefficients by FNV-64 mod4. Original active receipts point to the first unprocessed DFS node; all earlier subtrees in that shard are complete. For an unfinished shard, refine its remaining top-six prefixes by (((hash>>32)^hash)>>2) mod nsub. Only the owner reaches the saved partial subtree; the rest continue from its following sibling ancestry, respecting their secondary residue. The two original B shards and all 18 prescribed refinements terminated with EXHAUSTED, leaving no blocker; the independent receipt audit verified the whole partition. GSO lower bounds use downward rounded squared lengths and coefficients with propagated absolute integer error; actual leaves use integral coordinates/norm. T-order23 orbit reduction is valid on both sides: it permutes coordinates, preserves the affine class, and sends complete complementary fibres to complete fibres. No claim of half-depth from incomplete coverage.

**12-color tail.** Given a certified h, let H={0,h,Jh,J²h} modulo L. It is an order-four subgroup (J²+J+I=0, nonzero mod2 with no J-fixed nonzero vectors). All nonzero members have depth²≥21/2. Use CFR to choose b outside the four torus balls of squared radius 106743/10000 centered on H. Horizontal offsets H+{0,b,−Jb}. Differences between ternary groups are ±one of b,Jb,J²b plus H, each of depth²>106743/10000. Index-2 D4 parity-even site lattice on integer (a,b,c,d), color [(a+b+c) mod3,a mod2,(b+(a+b+c+d)/2) mod2], onto Z3×Z2². Its tail Gram is 65/24 on the first three diagonal, -1/24 between them, 21/8 in the last diagonal, zero others. For vertical integer squared length≤6, tail.py verifies every threshold exactly. Larger ordinary lengths have vertical squared length>16 because min eigenvalue 21/8 and ordinary square≥7. Zero vertical differences use parent minimum16. Rank-4 period is kernel of all colors, index12; its twelve sites cancel this index per center. Squared vertical covolume per center D=53361/256; prior exact D=10950866250772257481/50625000000000000. Radius2 occupied density≥6.367774598225213e−22, +1.870801239057239% from exact D ratio, periodic only. Parent theorem-backed existential offset b, not claimed explicit numeric coordinate or improved lattice lower bound.

## Replay

The [portable proof archive](proof-replay.tar.gz) contains the full source, pinned exact model and inputs, all 20 exhaustion receipts, and separate arithmetic/model checks. Archive SHA-256: `2d16ac908a226e6e3a959ff4bb4ffeec4b8e9891a8a365f73066993a4ff9175d`. A full independent rerun of enumeration is not claimed: an independent audit rederived the CRT/GSO inputs, verified the arithmetic ranges and complete receipt partition, and checked an interrupted mini-replay against a single uninterrupted enumeration.

```sh
mkdir -p replay && tar -xzf proof-replay.tar.gz -C replay
cd replay
(cd sphere36_mw92 && sha256sum -c FINAL_SHA256SUMS)
python3 sphere36_mw92/receipts.py sphere35_mw92 sphere36_mw92
python3 sphere36_mw92/audit.py sphere35_mw92
python3 sphere36_mw92/ranges.py sphere35_mw92
(cd sphere35_mw92 && python3 embed.py && python3 tail.py)
```

To recompute all eight original shards from scratch rather than trusting the saved checkpoint history, using a separate fresh output directory:

```sh
# Run from replay/ (GNU C++17, Python 3 with NumPy/SymPy for the audits)
g++ -O3 -std=c++17 sphere36_mw92/cover.cpp -o cover
mkdir -p fresh
for side in A B; do
  for shard in 0 1 2 3; do
    ./cover sphere35_mw92/$side.in sphere35_mw92/${side}d.in "$side" "$shard" 4 39 "fresh/$side.$shard.cover" 0 \
      >"fresh/$side.$shard.out" 2>"fresh/$side.$shard.log" || exit 1
    test ! -s "fresh/$side.$shard.out" || exit 1
    grep -q 'EXHAUSTED small minimum' "fresh/$side.$shard.log" || exit 1
  done
done
```

The checked floor in [result.json](result.json) uses rational bounds on π. The floor preserves the exact lower-bound direction despite decimal rounding.
