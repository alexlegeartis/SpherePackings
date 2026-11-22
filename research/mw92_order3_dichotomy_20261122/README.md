# MW92 order-three dichotomy: lattice lower bounds in 94–96D

| Dimension | Occupied density at least | Over previous repository lattice | Over previous unrestricted bound |
|---:|---:|---:|---:|
| 94 | 8.03567126556234e-21 | +3.2910% | +2.4390% |
| 95 | 2.10405317054093e-21 | +2.0335% | +1.7580% |
| 96 | 5.63315139086369e-22 | +2.4167% | No improvement |

Each new bound is attained or exceeded by a **single lattice**, not an
unrelated union of translates. These are existence proofs: they specify the
vertical forms in both cases and prove the existence of suitable horizontal
shears. They do not decide which case happens, or give explicit horizontal
offsets. The best previously certified **unrestricted** 96D bound remains
6.25082581067615e-22. The exact comparisons are in [results.json](results.json);
the named external baselines are the derived MW92 × A2 and MW92 × A3 products
and the 96D lattice in [Bos–Conway–Sloane (1982), Table 1](https://neilsloane.com/doc/Me103.pdf).

## Parent and torus depth

Take Shioda's even integral MW92 lattice L of minimum squared norm 16 and
covolume 47^7. Let J be its orthogonal order-three automorphism satisfying
I + J + J² = 0, as justified from the Mordell–Weil model and its height pairing
in the [previous construction](../mw92_extensions_20261122/README.md#parent-and-covering-input).
Put T = I − J. Write d(x) for the Euclidean distance from x to L. All the
covering estimates below are for **squared** distances.

The [Coxeter–Few–Rogers covering bound](https://doi.org/10.1112/S0025579300002059)
and the [cubic estimate](../adaptive_layers_20261115/COVERING_PROOF.md) imply
that any k prescribed centers on the torus leave a point at squared distance
strictly exceeding d_k from all centers. The exact rational d_k, for the values
of k used here, are in results.json. This follows since k Euclidean balls of
squared radius d_k per fundamental cell would have covering density smaller
than the universal lower bound. The same argument works on a quotient by an
index-three superlattice using a **single** ball with squared radius d_3:
lifting back to L would give three balls per L cell.

## The dichotomy

Select a with d(a)² > A = d_1. Put q = 16 − A and t = 16 − 3q.

- **If d(Ta)² ≥ t:** start with the vertical A2 Gram with diagonal q and
  off-diagonal q/2. Use horizontal shears a and −Ja. The first vertical shell,
  of squared length q, has horizontal differences in the orbit of a under
  ±J. The second shell, of squared length 3q, has differences in the orbit of
  Ta. All other vertical squared lengths are at least 4q > 16. These two
  assignments therefore suffice for the first two dimensions.
- **If d(Ta)² < t:** take a shortest representative w of Ta modulo L and
  write h = T⁻¹w and a0 = a − h. Then Ta0 belongs to L, and T scales lengths by
  √3 since TᵀT = 3I. In particular ||h||² < t/3. For a closest lattice point
  v to a0, apply v ↦ Jv + Ta0. This sends L onto itself, and the three relative
  vectors from a0 to the resulting lattice points have the same length and
  sum to zero. Average their squared distances from a0+h to obtain
  d(a)² ≤ d(a0)² + ||h||². Thus d(a0)² > A − t/3 = 32/3. For every u in L,
  the squared norm of a0−u is ||T(a0−u)||²/3, an even integer divided by 3.
  Hence d(a0)² ≥ 34/3. Moreover 3a0 belongs to L by (I−J²)T = 3I, and a0 is
  nonzero modulo L by its positive depth. Use a0 as the first shear with
  vertical squared length 14/3. For the second use a shear at squared
  distance > d_3 from L+Za0 (the index-three superlattice), vertical diagonal
  16−d_3, and vertical off-diagonal (14/3)/2. Every short difference with
  nonzero second coefficient has second coefficient ±1, so the required
  horizontal gap follows from this one choice on the superlattice quotient.

## Extend each branch and compute the bound

In 95D append one vertical basis vector from the larger Grams in
results.json, using d_12 in branch A or d_8 in branch B. In 96D, branch B
appends one more basis vector, using d_17. Orient each short vertical
difference with its last nonzero coefficient positive. The verifier exhausts
these vectors using inverse-Gram bounds, checks that this last coefficient is
1 at each extension step, and counts the forbidden centers. At the last
branch-B step, two coefficient lists differing only by a multiple of three
in the first coordinate give the same center. Each new shear can therefore
be chosen to miss all those centers to the corresponding depth.

### 96D branch A: Eisenstein multiplication

Use the 96D branch-A Gram in results.json instead of successive one-axis
extensions. Identify its two coordinate pairs with x,y in Z[J], and use
horizontal shift xa+yb. The squared vertical norm is

    5 (x0²−x0 x1+x1²) + (119/20) (y0²−y0 y1+y1²)
      + (4999/1000) (x0 y0+x1 y1−x1 y0).

For y=0, this has exactly the same two possible short shells as the first
branch, of squared lengths 5 and 15. Our a has squared depth exceeding A>11
and Ta has squared depth at least t>1, so these shells are safe. For y≠0,
multiplication by y on the torus is a surjective endomorphism: on real space
its squared length scale is y0²−y0 y1+y1²>0. Thus for uniform b, each event
d(xa+yb)²<16−Q(x,y) has measure at most the corresponding Euclidean ball
volume divided by 47^7, even if the ball overlaps its lattice translates.
Multiplying both x and y by the same unit in {±1, ±J, ±J²} preserves Q and
defines the same event. There are 13 short nonzero-y orbits, and the exact
sum of their ball-volume bounds is <0.997873010393<1. Hence there exists a
common b avoiding all the events.

All periodic horizontal translates are included by measuring torus distance.
For zero vertical difference, the minimum follows from L.

For each dimension the resulting lattice has squared minimum 16, attained
already in L, and squared covolume 47^14 times the determinant of the branch
Gram. Taking the **larger** of the two branch determinants gives the
unconditional bound. Its squared center density is 4^n divided by that
squared covolume; its occupied density is the center density multiplied by
the unit-ball volume. The verifier checks the displayed downward rounding
using a rational interval for pi.

From the repository root run:

```
python research/mw92_order3_dichotomy_20261122/verify.py
```
