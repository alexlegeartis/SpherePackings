# Exact limit for two arbitrary translated binary stocks in 13D

**Result.** Fix the cubic period `2 Z^13` at sphere radius 1. Take any
two binary length-13 codebooks of minimum Hamming distance 4, independently
translated by arbitrary real vectors in the torus. If their combined
occupancy is greater than 288, their cross squared separation is at most
`781/196 < 4`. Consequently at most 288 centers are possible in this
two-stock model, even with unequal occupancies and unrelated codes. The
sharp center density is `288/8192 = 9/256`, not a new overall record. The
classical `P13a` realizes it by the Julin 144-word stock at two heights,
with the horizontal all-half glide. [Leech–Sloane, §§2.6 and 4.3](https://doi.org/10.4153/CJM-1971-081-3);
[external density](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html);
[live 13D lower-bound table](https://github.com/alexlegeartis/SpherePackings/blob/main/LOWER_BOUNDS_1_96.md).

This does **not** bound three or more independently placed stocks, rotated
relative cubic periods, or variable point placements within one stock.
In particular it is not a proof of the optimum packing density in 13D.

## 1. A five-bit prefix lemma

Let A,B be binary length-13 distance-4 codebooks, with `|A|+|B|>288`.
On any five specified coordinates, there is a cross pair with Hamming
discrepancy at most one. This is true after independently XOR-translating
either full codebook, so it applies to any binary offset between stocks.

Suppose otherwise. Let S_A,S_B be their sets of projected five-bit
prefixes. They are disjoint and there is no hypercube edge between them.
Every prefix can carry at most `A(8,4)=16` suffixes. If either bank uses
at most two prefixes, its size is at most 32, while the other is at most
`A(13,4)=256`: contradiction. Otherwise take the smaller prefix set, of
size between 3 and 16. [Harper's vertex-isoperimetric theorem](https://doi.org/10.1016/S0021-9800(66)80059-5)
says its vertex boundary is at least that of an initial segment in weight
then lexicographic order. In Q5 these boundaries have sizes
`10,11,11,10,12,13,13,12,13,13,12,12,11,10` for sizes 3 through 16.
Since the other bank avoids this entire boundary, the union S has at most
22 prefixes. The checker independently obtains the boundary list without
searching arbitrary cube subsets.

Take a maximum matching in the subgraph of Q5 induced by S. Every edge
lies within one bank. An adjacent pair has combined occupancy at most
`A(9,4)=20`, by keeping the differing prefix coordinate with the eight
suffix coordinates. An independent set of Q5 has size at most 16, by a
coordinate perfect matching. By König's theorem the induced bipartite
matching therefore has size at least `|S|-16`. The total occupancy is
at most `16|S|-12(|S|-16) = 4|S|+192 <=280`, again impossible.

The binary-code cardinalities used throughout are the exact entries in
[Brouwer's coding table](https://aeb.win.tue.nl/codes/binary-1.html),
including its original references. `A(12,3)=A(13,4)=256` follows from
parity extension and puncturing. Standalone upper inputs are `A(8,4)=16`,
`A(9,4)=20`, `A(11,4)=72`, and `A(13,4)=256`.

## 2. From three combinatorial constraints to all real glides

In the cubic period choose a nearest binary corner to the real glide
coordinatewise, absorb it in the second stock, and sort the coordinates
by the mismatch penalty `0<=a1<=...<=a13<=1`. If the real coordinate is
at distance `(1-ai)/2` from its nearest integer, the contribution to
squared separation is `(1-ai)^2/4` for a matching bit and adds `ai` for
a mismatching bit. Integer torus reduction is exactly coordinatewise.
Let `f(a)` be the minimum total mismatch penalty over cross codeword pairs.
If total occupancy exceeds 288, the following three bounds hold:

1. **Puncture the lightest coordinate.** Each punctured code has minimum
   distance 3 and unchanged size. If all cross projected discrepancies
   were at least 3, their disjoint union would exceed `A(12,3)=256`.
   Therefore `f <= a1+a12+a13`.
2. **Fix the two heaviest coordinates.** Some prefix cell has combined
   occupancy greater than 72, since there are four cells. Each bank has
   at most `A(11,4)=72` suffixes in that cell, so both occur. The cross
   suffixes cannot all be bitwise complements: then each nonempty set
   would have size one. Thus a cross pair also agrees in at least one of
   the remaining coordinates, giving `f <= sum(a2,...,a11)`.
3. **Use the five-bit prefix lemma on the five heaviest coordinates.**
   It gives `f <= sum(a1,...,a8)+a13`.

Adding the common base `sum_i (1-ai)^2/4` to the minimum of these
three linear bounds gives a maximum of exactly `781/196`. For
completeness, the domain `0<=a1<=...<=an<=1` is the simplex with vertices
`v_k` consisting of k ones at the end. Split it according to which bound
is smallest. On each resulting polytope the upper objective is convex,
so its maximum is attained at a vertex. These vertices are a simplex
vertex, an edge cut by one equality, or a triangle face cut by the two
independent equalities. Degenerate equalities introduce no extra case:
a vertex with s positive simplex coefficients needs at least s-1
independent active equalities. The exact rational check `verify.py`
includes all these candidates, even those not satisfying a particular
region's inequalities, and evaluates the lower envelope at each. Its
maximum is the vertex `(v13+5 v3+v2)/7`.

The bound `781/196` is on the particular upper envelope; we do not claim
it is an attained cross separation of a literal pair of codebooks.

## 3. Independent 12D heterogeneous whole-row consequence

If two binary length-12 distance-4 codebooks have combined occupancy
greater than 256, their cross squared separation under any real cubic
glide is at most 3. This holds for **unrelated** and **unequal** stocks,
not only the fixed 27 unit-row banks of the earlier 13D work.

First, after any binary shift, the cross mismatch count is at most 2;
otherwise adjoining opposite binary marker bits gives a length-13
distance-4 code exceeding 256. Secondly, any prescribed seven coordinates
have a matching cross pair. Each projection cell has at most
`A(5,4)=2` suffixes, so the sum of projected cardinalities is larger than
128, forcing intersection. For `a1<=...<=a12` these give respectively
`a11+a12` and `sum(a1,...,a5)` as upper bounds for f. The same exact
simplex argument, now with only one cutting hyperplane, gives maximum
3, attained only at `a=0` (the all-half glide). In 13D this forces a
vertical gap of at least 1 between any such pair of whole 12D rows.
It does not claim the same gap for unequal rows whose combined
occupancy is at most 256.

## A constructive lead and its scope

The main theorem leaves three-bank placement or a changed period metric
as precise escape routes. A first candidate palette chooses three
half-grid offsets on the last 12 coordinates of 13D, grouped into three
four-coordinate blocks: `0`, `1/2` on blocks 1+2, and `1/2` on blocks 1+3.
For each pair of banks the horizontal half-grid contribution is exactly 2;
the remaining five bits must differ in at least two places. The three pairs
use three distinct five-bit projections. This is a different problem from
the boxed two-bank proof, and no improving motif for it is claimed here.

Run `python3 verify.py` for the exact simplex and Harper arithmetic. No
catalogue of actual codewords or numerical geometric search is used by
the proofs.
