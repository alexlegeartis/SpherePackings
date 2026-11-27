# The standard 13D Julin packing is saturated even against arbitrary real points

This is a point-insertion barrier for the particular P13a stock used by
the [repository's rigidity proof](https://github.com/alexlegeartis/SpherePackings/tree/main/research/julin13_rigidity_20261126),
not a prohibition on deleting or rearranging any old centers. It uses
neither the two-stock theorem nor a numerical hole grid.

Take the explicit 144-word Julin code J on twelve bits from
[Leech–Sloane, §2.6](https://doi.org/10.4153/CJM-1971-081-3), with the
six loose pairs at `(0,1),(2,3),...,(10,11)` as in that repository proof.
In the horizontal cubic torus modulo 2, its squared covering radius is
at most `103/32`. Indeed every projection onto any retained set of
`5, 9, 11, 12` coordinates has Hamming covering radius at most
`0, 1, 2, 3`, respectively. `verify_saturation.py` constructs the code
from the eleven cyclic shifts, and checks these four facts independently
by XOR-expanding the projected words with all masks up to the stated
radius. They are coverage certificates, not numerical approximate CVPs.

Given a real point, fold to a nearest binary corner and sort mismatch
weights `a1<=...<=a12` as in the main note. A projection covering bound
`r` on the `k` heaviest coordinates provides a word with at most `r`
mismatches there, whatever its remaining coordinates. Thus the weighted
mismatch penalty is bounded by the minimum of the following four:

```
sum(a1,...,a7),
sum(a1,a2,a3)+a12,
a1+a11+a12,
a10+a11+a12.
```

Adding the common base `sum_i (1-ai)^2/4`, the same simplex argument as
in the main note applies, now cutting with up to three active equalities.
The checker solves the resulting systems over exact rationals. The
largest upper-envelope value is `103/32` (e.g. at ten 1/4 weights,
one 1/2, and one 1). No equality for the literal covering radius is claimed.

The 13D packing has a full translated copy of this stock at every
integer height, with glide alternating between 0 and the all-half
horizontal vector. For an arbitrary real point choose the nearest
integer height, at distance at most 1/2, and then the closest horizontal
stock center there. The squared distance is at most `103/32+1/4=111/32<4`.
Every new sphere of radius 1 therefore overlaps an existing one. Any
improvement of this particular stock must remove or move old centers,
not merely insert in a previously unused glide phase.

Replay: `python3 verify_saturation.py`.
