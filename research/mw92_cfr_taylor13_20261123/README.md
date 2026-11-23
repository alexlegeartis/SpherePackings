# Exact higher-order covering refinements for MW92

This package improves prior certified bounds strictly but only by very small amounts:

| Dimension | Packing class | Occupied density at least | Increase over prior repository guarantee |
|---:|:---|---:|---:|
| 93 | Lattice and unrestricted | 3.045564894259504e-20 | 0.00004977% |
| 94 | Lattice and unrestricted | 8.036027386567117e-21 | 0.00443175% |
| 95 | Lattice and unrestricted | 2.104056220143222e-21 | 0.00014494% |
| 96 | Unrestricted periodic only | 6.250833919802144e-22 | 0.00012973% |

The 96D lattice result does not change. All horizontal offsets are existential.

## Covering lemma

Use the same even MW92 parent of minimum squared norm 16 and covolume 47^7, and its integral orthogonal J satisfying J²+J+I=0, as in [the parent construction](../mw92_extensions_20261122/README.md). Let U be uniform on the simplex with 92 coordinates, S=sum U_i² and mu=2/93. The [previous covering proof](../adaptive_layers_20261115/COVERING_PROOF.md) gives a Coxeter–Few–Rogers lower bound 93(92/93)^46 E[(1+S)^−46]. The degree-13 Taylor polynomial of (1+s)^−46 about mu is a pointwise lower bound on s≥0: its degree-14 Lagrange remainder is nonnegative even for s<mu. Using Dirichlet integration and the multinomial formula gives exactly

    E[S^j] = j! [X^j] (sum_{i=0}^j ((2i)!/i!) X^i)^92 / (92)_(2j).

Here (92)_(2j) is a rising factorial. `moments.py` computes raw and centered moments and evaluates this rational lower bound. It strictly improves the analogous degree-3 value. The verifier proves that any union of k torus balls of squared radius d[k] fails to cover by exact comparison of k·pi^46·d[k]^46/(46!·47^7) with the new bound, taking the upper endpoint for pi. For k=1,3,7,8,12, the respective d[k] are 11.00098526, 10.74136270, 10.54532299, 10.51475582, 10.42248121.

## Vertical forms and exact comparator

Set p=16−d[1], t=16−3p. Squared minimum 16 gives squared center density 4^n/(47^14 D), where D is the squared vertical covolume per center. The strict ordering is certified by comparing exact rational D, not rounded occupied densities. These are the old and new values:

| n | Old D | New D |
|---:|---:|---:|
| 93 | 1249754929/250000000 | 249950737/50000000 |
| 94 | 85939/4500 | 859313833/45000000 |
| 95 | 4581621816246316217474355911/62500000000000000000000000 | 73305736561980600794262547/1000000000000000000000000 |
| 96 | 2737723665915453601/12656250000000000 | 10950866250772257481/50625000000000000 |

* In 93D use vertical Gram [p] with one shear a at squared torus depth greater than d[1].
* In 94D and 95D use [the existing order-three dichotomy](../mw92_order3_dichotomy_20261122/README.md) with the new depths. With T=I−J, if d(Ta)^2≥t, the first two shears are a, −Ja. The first two diagonals are p; the third diagonal, if used, is 16−d[12]. If d(Ta)^2<t, apply the old averaging argument to get an exact order-three fixed class a0 with d(a0)^2≥34/3 (the cutoff d[1]−t/3=32/3 is unchanged). The diagonals are 14/3, 16−d[3], 16−d[8], truncated if appropriate. In both cases off-diagonal i<j is half diagonal i. The rank-2/3 short-vector enumeration checks all remaining torus choices and the maximum over the two branch determinants yields the guaranteed density.
* In 96D use vertical parity-even integer D4 coordinates w=(w0,w1,w2,w3), metric z=8/3, s=(3(16−d[7])−2z)/4, and norm
  `s w0² + z(w1²+w2²+w3²) + (s−z)(w1+w2+w3)²/3`.
  The four eigenvalues in the ambient integer coordinates are s,s,z,z and the underlying D4 has index two, so D=4s²z². Use colors `(w0+w2+2w3, w1+w2+w3) mod 3`. This map is onto; the nine translates of its kernel cancel index nine in the period covolume. Choose horizontal torus offsets by [the inherited nine-offset lemma](../mw92_extensions_20261122/README.md#nine-translates-in-96d). If the second difference color is nonzero then horizontal squared separation exceeds d[7], and if only the first is nonzero then it exceeds d[1]. The checker excludes a dangerous zero-color difference or any insufficiently long nonzero-color difference. Since the least eigenvalue is z, checking integer vectors of ordinary squared length below seven is sufficient. The zero vertical difference is safe in all constructions by the parent minimum.

## Exact replay and pins

Run `python -B verify.py` here; the expected output is `result.json`. The
checker is standard-library-only and requires Python without `-O`. The
two copied exact previous result files are integral inputs. Both bytes were checked
against the repository Git blobs:

* `old.json` = `research/mw92_extensions_20261122/results.json`, Git blob `a254e4c8494a20c9feb74f6e2edf45fb09f89843`.
* `prior95.json` = `research/mw92_order3_dichotomy_20261122/results.json`, Git blob `b8834f3442381576471da20062188595914eced8`.

The separate `SHA256SUMS` pins the package.
