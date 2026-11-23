# A 48D obstruction for the zero-sum E8 tensor carrier

Let `L = {(u0,…,u6) in E8^7 : sum(ui)=0}/sqrt(2)`, with the Euclidean
metric inherited from the sum-zero subspace. Any union of translates of
`2L` with representatives in `L` has center density at most
**38,246,144 / 2401 ≈ 15,929.2562**, below the P48 comparator
`(3/2)^24 ≈ 16,834.1122`. Neither additivity nor fullness of the chosen
motif is required. This is an upper bound for this carrier, not for
arbitrary 48D packings, and it does not claim attainment.

## Proof

The covolume of the unscaled `A6 ⊗ E8` is `7^4), so the period `2L`
has covolume `2^24·7^4`. The period contains a vector of squared norm 8
from two opposite roots, and every displacement in `L) has integral
squared norm. Thus an actual squared minimum `g` lies in {1,…,8}.
A motif of size N has center density `N·(g/8)^24 / 7^4`.
For g ≤ 3 the trivial N ≤ 2^48 suffices.

The labels of `L/2L` are seven elements of `V = E8/2E8` of sum zero.
This model is onto: any lift of the labels has sum in 2E8, which can be
removed from one lift without changing its label. Split by the first
label (256 choices). Within each fiber the remaining six-coordinate
difference labels sum to zero and form a group H of size 2^40.
The verifier checks the root dictionary in doubled E8 coordinates and
lifts the full diagonal O⁺(8,2) action using explicit reflection
stabilizer orbits of degrees 135, 64, 35, 16, 9, 4, 2, 1. Together
with position permutations these preserve the joint metric. Using
[Witt extension](https://web.math.ku.dk/~wahl/Forms.pdf), the ordered
orbits are given by a kernel and the quadratic form on its quotient.
An orbit is forbidden only when the integer enumeration exhibits an
actual representative with sum zero and squared norm below g. Leaving
an orbit allowed when the enumeration misses a short lift is safe.

The difference distribution of a fiber has value 1 at zero and
nonnegative Fourier transform, even after orbit averaging. If N_j is
the size of allowed nonzero orbit j and M_ij is the corresponding
character sum, its mass a_j satisfies sum_j(-M_ij/N_j)a_j ≤ 1
for each nontrivial character orbit i. The code generates the integer
sums with a Walsh transform over quotient quadratic forms and checks
their aggregate counts. The supplied nonnegative rational LP dual
checks every column and bounds the motif by 256 times the integer fiber
cap. All five dual checks and comparisons use exact arithmetic.

| Squared min. | Fiber cap | Full motif cap |
|---:|---:|---:|
| 4 | 4,295,571,857 | 1,099,666,395,392 |
| 5 | 380,673,372 | 97,452,383,232 |
| 6 | 31,318,314 | 8,017,488,384 |
| 7 | 2,031,616 | 520,093,696 |
| 8 | 149,399 | 38,246,144 |

The maximum center bound is achieved by the *certificate* at g=8.
The same certificate excludes a record from common linear
deformations with Gram condition number at most 20047/20000:
after rescaling the bottom eigenvalue to 1, if the top is K the
squared minimum increases by at most K and the Gram determinant is
at least K, so center density grows by at most K^(47/2).

## Replay

Extract [replay.tar.gz](replay.tar.gz) and run `python3 verify.py` with
GCC supporting C++20. It regenerates the five integer Walsh tables,
checks their SHA-256 hashes against the frozen exact replay, and
verifies the rational duals and final density comparisons. The
generated multi-megabyte tables are not required as inputs.

Archived replay SHA-256: `9aa0b7d2475c319b387a58bedcfee72d5dd93a368b3190910ff3a3a56a460251`.
