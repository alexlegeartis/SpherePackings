# E8 tensor center-density cap: two heavy slots

In the 48-dimensional seven-slot zero-sum E8 half-period model, a
positive-definite weight tuple `(a,b,c,1,1,1,1)` with `a<0` and
`b,c>=2` has center density at most **16,800**, for any half-period
motif, linear or nonlinear. This is below the center density
`(3/2)^24` of the 48D record lattice. It is not an unrestricted
packing bound.

Use the normalization and Fourier lift lemmas in the accompanying
[Round 20 determinant proof](../../scoped_density_round20_20261124/48-e8-determinant-negative). In particular, with
`D=ab+ac+bc+4abc`, squared minimum `g`, and label count `N`,
the density is `N(g/8)^24/D^4`. The determinant-only Hermite and
period bounds cover `D<=11/4000` and `D>=6`.

For fixed D the value `a=(D-bc)/(b+c+4bc)` decreases in either
heavy variable: its b derivative is
`-(c^2+D(1+4c))/(b+c+4bc)^2`, and its c derivative is symmetric.
Consequently upper bounds on literal lifts and periods involving
only a and ordinary weight-1 slots, computed at b=c=2,
remain valid on the entire unbounded domain b,c>=2. The
two-slot shortcut uses a+1, no greater than a+b. Shapes 31 and 41
change only a and ordinary slots. The finite proof enforces these
conditions and allows the uniform squared-minimum bounds 4 and
8 times the cube root of D. It is valid even for fictitious
reference endpoints where the recomputed a is positive: only
the upper bounds on literal lengths, not a feasible reference
packing, are used.

The twelve determinant slabs form a complete cover, and the
inner cut lists split the complete range from 0 to the supplied
squared-minimum bound. Fourier dual weights are nonnegative
integers and all inequalities are recomputed as exact rationals.
The discovery search does not participate in replay.

Run with a standard C++20 compiler and Python 3:

    python3 -B replay.py proof.json.gz

Combined with the earlier three-exception proofs for all
nonnegative weights, all exceptional weights at most 1, and
the negative chart `0<b<=2,c>=max(1,b)`, this closes the full
ordered three-exception family.
