# 48D E8 half-period motifs: negative slot and second exceptional weight at most two

Let the carrier be the sum-zero sublattice of seven copies of E8 scaled by 1/sqrt(2), with squared slot weights (a,b,c,1,1,1,1). Let a motif be an arbitrary union of cosets of twice this carrier, without additivity assumptions. **For a<0, 0<b≤2, c≥1, b≤c, and positive definite scalar metric, the center density is at most 16800.** There is no upper restriction on c. This closes a new negative/heavy region, including the near-cancellation boxes that were open in the original (a,b,c)-coordinate scan. It excludes improvement over the P48 density (3/2)^24 in this chart only. It is not an upper bound for unrestricted 48D packings, and does not address the intermediate case b>2 unless covered by other results.

## Determinant and length parametrization

The squared-slot scalar determinant on the six independent slots is

    D=ab+ac+bc+4abc .

Consequently a=(D−bc)/(b+c+4bc). The 48D period covolume is 2^24 D^4, and N labels at squared minimum g give center density N(g/8)^24 / D^4. Distinct labels number at most 2^48. The two-slot opposite-root period gives g≤4(1+a)<4. The three-dimensional period sublattice in a single E8 root direction with coefficients (−(u+v+4t),u,v,t,t,t,t) has scalar Gram determinant 4D before the period norm factor 4. The classical rank-three Hermite bound gamma_3=2^(1/3) therefore gives g≤8D^(1/3). In each determinant interval [D0,D1] it supplies the separable bound

    center density / N ≤ min[(g/8)^24/D0^4, (g/8)^12].

In particular D≤11/4000 is closed by 2^48 D^4<16800; D≥6 is closed by 2^24/6^4<16800. The package covers every determinant between those boundaries.

At fixed nonnegative D, the expression for a is decreasing and convex separately in b,c: for example a_b=−[c²+D(1+4c)]/(b+c+4bc)^2. It is affine in D. Every weighted squared length used below has a nonnegative coefficient of a, and hence its maximum over any box in (D,b,c) is attained at a corner (taking a single fixed physical witness). This is stronger near determinant cancellation than independently upper-bounding original a,b,c at an unrelated corner and estimating D at an opposite one.

## Period and fibre tests

A root-direction period vector is chosen with integer coefficient u at b, v at c and ordinary coefficients as equally distributed as possible with total t, then coefficient −(u+v+t) at a. Its squared length is exactly

    4[a(u+v+t)^2 + b u^2 + c v^2 + 4 floor(t/4)^2 + (2 floor(t/4)+1)(t mod 4)].

A saved vector norm evaluated by the exact corner maximum or the rank-three Hermite bound certifies the g upper bound for the tile. The verifier checks that every saved coefficient tuple belongs to the literal period set and checks the appropriate algebraic inequality for any rounded cube-root cap.

Two-slot phase fibres of weight-sum w=a+b or a+1 have size at most 16 for g>w and at most 1 for g>2w. The differences of two phases in such a fibre are labels in E8/2E8, with a literal opposite lift of carrier squared length w for norm-2 classes and 2w for nonzero isotropic classes. In the first range all pairwise differences are isotropic, and after translating to include zero the span is totally singular in dimension 8, hence rank at most 4. This is the source for the inexpensive pair clauses in the tree.

All other fibre clauses use nonnegative exact Fourier character certificates on bounded changing-slot groups, inherited from the previous zero-sum half-period implementation. For table weight code 22/31/32/40/41/42 choose the stated k+1 slots (first digit k), with respectively the first two exceptional weights, only a, both a,b, no exceptions, only a, or both a,b. The code 40 uses four ordinary slots and a fifth slot a≤0, treated as ordinary at upper weight 1. Other fibres use the literal four ordinary slots as needed. Freeze the other 6−k independent labels, so the global count bound is the fibre cap times 256^(6−k).

For each nonzero difference orbit choose any listed actual short lift. For each such lift bound its weighted squared length by the maximum over the eight box corners, first rounding the scalar coordinates upward to a rational grid of denominator 10^8. Minimize this upper bound over the listed lifts. If it is at most the lower endpoint l of a range l<g≤h, that difference orbit is forbidden throughout the range. This convention covers thresholds without any floating-point tolerance. Given exact orbit cardinalities n_j and character sums W_ij, the verifier checks the nonnegative integer multipliers v_i. If z=min_j[-sum_i v_i W_ij/n_j]>0 over remaining nonzero orbits, the fibre cap is floor[1+sum_i v_i/z] by summing the positive-definite character tests. The proof is valid for nonlinear motifs. A missing/no-sign certificate can only revert to the full finite-fibre cap, never infer exclusion.

## Why the heaviest slot may be arbitrarily heavy

The finite exact cover splits 11/4000≤D≤6, 0≤b≤2, 1≤c≤4. For each leaf touching c=4, the verifier checks that its period bound follows either from g≤4, from the Hermite ceiling with exact D upper endpoint, or from a literal period witness whose coefficient on the c slot is zero. With D,b fixed, a decreases as c increases, so that period witness only gets shorter. Similarly all fibre clauses touching c=4 freeze the heaviest label; their other weights stay unchanged except a, which decreases. The same fibre bound, replayed at c=4 exactly, remains valid for all c≥4. The one empty tile at c=4 reaching the tail is instead certified directly by g≤4b. This extrapolation makes no assertion about c-dependent witnesses in interior tiles. Together with the universal small/large D cuts, it closes all positive determinant values for the stated infinite tail.

## Exact replay

Run `python3 replay.py proof_4_2_tail.json.gz`. Python standard library, a C++20 compiler and a little-endian host are needed; there is no numerical optimizer, input downloads, or Git step. The independent replay rechecks the root dictionary and lifted orthogonal group, regenerates only the Fourier tables actually used, and verifies every exact leaf and every g subinterval. The initial determinant slabs form a dyadic partition, but internal tree splits are chosen adaptively. The compact box uses 404 bisections (367 fibre/paired proof leaves, 40 period-only, 9 empty); the tail replay also checks all 57 leaves touching c=4. No open leaf is accepted. Source root_dictionary.py is bundled for the pinned E8 root dictionary verification.
