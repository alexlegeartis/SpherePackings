# 48D seven-slot zero-sum E8 gluing: exclusion of every one-exceptional metric, including the dual tensor

**Result.** Give the seven slots of the zero-sum carrier
`L = { (u_0,...,u_6) in E8^7 : sum u_i = 0 } / sqrt(2)`
squared metric weights `(a,1,1,1,1,1,1)`. This is a positive-definite
metric on the zero-sum subspace precisely when the real number
`a > -1/6`, even though the diagonal weight in the ambient space
can be negative. An arbitrary union of cosets of `2L` with
representatives in `L` has center density at most
`38,246,144/2401 < 15,930 < (3/2)^24`. For `-1/6<a<=0` the tighter cap
is `1,213,025,622,610,333,925,376/79,792,266,297,612,001 < 15,203`.
Neither bound claims attainment.
The scope includes every nonlinear or additive half-period motif
in this full family, in particular on the dual tensor
`(A6* tensor E8)/sqrt(2)` at `a=-1/7`. It does not cover other
periods, more independent weight ratios, off-diagonal metrics,
or other 48D packings.

The comparator `(3/2)^24` is attained by the extremal 48D even unimodular lattices listed by [Cohn](https://cohn.mit.edu/sphere-packing/) and in [Gonçalves–Vedana, Theorem 1](https://doi.org/10.1017/fms.2025.9). Their exact global assertion has a specific forbidden-distance hypothesis; no such global claim is made here.

## Geometry and density

Normalize the six common weights to 1 by an overall dilation. Eliminating the exceptional slot leaves scalar Gram `I_6+aJ_6`, so the covolume of `2L` is `2^24 (1+6a)^4`. There is a literal period contact of square `4(1+a)` using opposite roots in the light and an ordinary slot, one of square 8 using two ordinary slots, and one of square `24(1+6a)` by repeating the same root in all six ordinary coordinates and using its minus-sixfold multiple in the exceptional slot. Let `g` be the actual squared minimum and `N` the number of distinct phases. The center density is `N (g/8)^24 / (1+6a)^4`.

Identify `V=E8/2E8` with the binary quadratic space of dimension 8, with q determined by half the original squared norm modulo 2. A sum-zero label tuple lifts to a sum-zero E8 tuple: subtract the sum of any lifts from one coordinate. The tuple labels therefore identify `L/2L` bijectively. Fixing three ordinary labels partitions the motif into at most `2^24` fibres; differences then belong to a group of four zero-sum slots (one exceptional and three ordinary), with three free labels. We always fix ordinary labels for `a<=1`; in the negative chart we also use three-slot fibres by fixing four ordinary labels (at most `2^32` fibres). At `a=1` all slots are interchangeable. For `a>=1`, fix the exceptional label and two ordinary labels instead, so all four changing labels are ordinary and their metric is exactly the one at `a=1`. The separate six-slot bound below also uses only ordinary slots. This proves the same motif-size caps as at `a=1` for all `a>=1`; the density denominator only increases.

At `a=-1/7` the scalar Gram `I_6-J_6/7` is the inverse Gram
of the ordinary `A6` basis `e_i-e_0`, so this is exactly the dual
tensor with the same normalization (E8 is self-dual).

## Weighted three- and four-slot witnesses

The solver-free replay validates a root-class dictionary and generates the full diagonal orthogonal-group orbits by a kernel and a pulled-back quadratic form (Witt extension in the split nonsingular quadratic space; see [Witt's lemma](https://web.math.ku.dk/~wahl/Forms.pdf)). Reflection-generated stabilizer orbits verify all the label isometries used have E8 lifts. Permuting the three ordinary slots is also an isometry for all a.

On each orbit we find literal lifts of score `h+a*l` with integer nonnegative scores `h,l`. Choose a forced slot and shortest E8 representatives in the other occupied slots, using the zero vector for their zero classes. The forced vector, which may be longer than its own class minimum, is minus the sum. In doubled E8 coordinates a length score is a squared integer norm divided by 8. A witness is sufficient to forbid its whole orbit when its score is strictly less than `g`; completeness of this witness search is unnecessary. For each heavy score the smallest light score is retained on the
nonnegative chart. Pruning a witness when it is never strictly below
the period ceiling on `[0,1]` is optional for validity. On the negative
chart always force the exceptional (light) coordinate and choose the
largest light score instead, because that makes the total score smaller.
Every retained witness is a literal lift in the positive-definite
zero-sum metric.

The `10 x 10` and `44 x 44` integer Fourier matrices count characters on the nonzero difference orbits. One can first count with a single pulled-back form without fixing the exact kernel. The inverse Walsh transform gives these counts: the sum for one hyperbolic coordinate pair, with k-bit coefficients and form mask M, is
`sum_{x,y} (-1)^(x^T M y + A*y + B*x) = 2^k sum_{x:M x=A} (-1)^(B*x)`.
Multiply the four hyperbolic-coordinate sums. Subtract counts of strictly larger kernels, then fuse ordinary-position permutation orbits. The dual uses the two or three independent character labels after the exceptional one is set to zero by diagonal translation. A character of a code difference distribution has nonnegative Fourier coefficient, being a squared modulus after normalization. If orbit j has cardinality `N_j` and `sum_char(i,j)=W_ij`, any nonnegative row weights `v_i` with
`z = min_j(-sum_i v_i W_ij / N_j) > 0` over all remaining nonzero orbits
prove that the fibre size is at most `floor(1 + sum_i v_i / z)`.
The replay verifies every such inequality using exact fractions.

For the high window `g>7` there is also the earlier six-ordinary-slot Fourier proof, fixing just the exceptional label. Uniform short lifts of squared norm at most 7 suffice. The replay independently regenerates the corresponding matrix and checks the archived exact rational dual: the fibre bound is 149,399, hence `N <= 38,246,144`. This part uses the same Fourier argument and does not need a second solver.

## Finite coverage of the continuum

Between intersections of the recorded affine witness lines, integers 0 through 8, and the period lines `4+4a` and (on the negative chart) `24(1+6a)`, their order is constant. For adjacent lines L and U, take `L(a) < g <= U(a)`. A witness at or below L at a rational midpoint remains at or below L throughout the closed parameter interval, including its endpoints. Thus the disallowed orbits, and their verified bound, apply across the whole window. Equality to a lower line is covered by the preceding window. Degenerate windows at intersections are empty. For the nonnegative chart there are six parameter intervals and 96 windows; for the negative chart there are 45 intervals and 701 windows. Trivial phase cardinality and, for the nonnegative chart, the separate uniform six-slot dual are also used. For any affine upper line `U(a)=h+l*a`, the derivative sign of `U(a)^24/(1+6a)^4` wherever positive is that of `l-h+5al`. This is nondecreasing, so the maximum is at an endpoint. The verifier uses both endpoints exactly without discretizing any parameter or relying on approximate solver status.

At the excluded endpoint `a=-1/6`, the density factor is
interpreted by its right limit zero: every present window has
`0<=U(a)<=24(1+6a)`. The replay checks the stated exact smaller cap on the negative chart.

Replay: `python3 verify.py` in this directory (Python stdlib and g++ C++20;
no solver or network). The binary matrix format used only inside the
temporary replay is little-endian.
