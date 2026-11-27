# Four optical colors over FCC in dimension 15
Date: 2026-11-27. No improved density is claimed. The exact
[15D laminated benchmark](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html)
is \(1/(16\sqrt2)\) in center-density units.

## Four-color construction and sharp unit-metric cap

Use the FCC normal set \(D_3/\sqrt2\), with normal points indexed
by \(x\in\mathbb Z^3\) of even coordinate sum. Let \(c_i=x_i\bmod2\).
Split twelve horizontal binary coordinates into blocks
\(B_0,B_1,B_2\) of length four. The horizontal period is \(2\mathbb Z^{12}\),
and the optical offset at \(x\) has value \(c_i/2\) in block \(B_i\).
There are four offset colors, rather than the binary pair of offsets
zero and the all-half vector. At a normal edge of unit length, two
normal coordinate parities differ; the optical offsets differ by
one-half in those eight horizontal coordinates. These coordinates
contribute exactly 2 to horizontal torus squared separation, regardless
of the binary words. In the remaining common block, at least one bit
must differ: the two projected row codes on that block must be disjoint.

A row motif of minimum squared separation 4 is a binary distance-four
code of length twelve. The upper bound 144 is the consequence of
\(A_2(10,3)=72\) in the
[Östergård–Baicheva–Kolev proof data](https://users.aalto.fi/~pat/72.html),
by shortening once and then puncturing once. If its three projection
sizes are \(a,b,c\), fixing two block patterns leaves at most two
words (only an antipodal pair can have Hamming distance four in the
remaining block). So its size is at most
\(\min(144,2\min(ab,ac,bc))\).

Every unit FCC tetrahedron contains all four colors. Its opposite
edge pairs give the three different perfect matchings on the
four vertices; each matching constrains disjoint projection sizes
in one specified block. For each matching the two pairs of projection
sizes sum to at most 16. Enlarge each pair to sum exactly 16 without
decreasing the bound. The exhaustive integer comparison in
verify.js evaluates the \(17^6\) possibilities for the six sums
and gives a cap of 512 total row words per tetrahedron.
For example orient the fixed tetrahedron with raw positions
\((0,0,0),(0,1,1),(1,0,1),(1,1,0)\). Averaging its translates
gives a cap of 128 words per FCC site, even with varying rows
and no period. This gives the density \(1/(16\sqrt2)\), not a gain.
This first argument assumes the four indicated optical masks and
unit FCC edge lengths.

Here is exact attainment in this model. Number coordinates in
each block by the two bits \(b\in\mathbb F_2^2\). Let \(C\) be
the kernel of the five-bit parity check with columns \((e_i,b)\)
on block \(i\). The columns are distinct and nonzero; no three
sum to zero by the block parities. They span rank five, so
\(C\) is a binary \([12,7,4]\) code of size 128.
At \(x\), use its coset whose block parities are
\(p_i=c_{i+1\bmod3}\), and whose last two syndrome bits encode
\(\sum_i\lfloor x_i/2\rfloor\bmod2\) as a fixed nonzero-bit toggle.
At a unit edge the common block parities differ. At an axial
two-step, normal squared separation is 2, the optical and three
block parities are unchanged, and the toggle gives binary distance
at least 2. At a normal squared separation of 3, the masks differ
on eight coordinates, more than enough. Other normal differences
clear the packing by normal distance alone. The exact periodic
check covers all literal differences with normal squared separation
below 4.

## Real covering radius of the whole horizontal layer

Let \(H=C+2\mathbb Z^{12}\), of covolume 32, so one complete
horizontal row is a translate of \(H\). Its real covering
radius squared is exactly 3, not merely a half-grid bound.
To prove this, round a real point to its nearest integer coordinates,
and write \(r_j\in[0,1/2]\) for the magnitudes of its residuals.
Changing a rounded bit adds a cost \(w_j=1-2r_j\) to the starting
squared distance \(R=\sum r_j^2\). To obtain a point of \(H\),
the changes must have the syndrome of the rounded word.

* If two or three block parities are nonzero, correct by one
  change in every such block, uniformly over choices with the
  required summed two-bit syndrome. Every change is uniform
  separately. The expected total in a used block is
  \(\sum r_j^2+\tfrac14\sum w_j\le1\), and the cost in another
  block is at most 1.
* If the block parities vanish and the two-bit syndrome is nonzero,
  it can be corrected by one of six bit pairs (two disjoint pairs
  in each block). Averaging gives total at most
  \(2+\sum_j(r_j^2-r_j/3)\le3\).
* If just one block parity is nonzero, change a bit of largest
  residual \(u\) in that block. Its cost including rounding is
  at most \(4u^2+1-2u\le1\). If needed, correct the remaining
  two-bit syndrome in one of the two inactive blocks using the
  cheaper of two complementary pairs there. With total change
  weights \(W_1,W_2\) in those blocks, the sum of their rounded
  costs is at most \(2-(W_1+W_2)/4\), and the extra pair cost
  is at most \(\min(W_1,W_2)/2\). Thus the total is at most 3.

The zero syndrome needs no change, and the integer syndrome
with all three block parities nonzero attains squared depth 3.

## Coupled metric: exact aperture and an octahedron obstruction

Allow one full translated copy of \(\sqrt t\,H\) above each
point of \(\sqrt{\gamma/2}\,D_3\), at unit sphere radius.
Here \(\gamma\) is the normal FCC edge squared length, and \(t\ge1\)
is the multiplier of horizontal squared lengths. There is no
linearity or periodicity condition on the shift assignments.
The center density divided by the benchmark is
\(t^{-6}\gamma^{-3/2}\), so a proposed gain must satisfy
\(\gamma<t^{-4}\).

A necessary horizontal torus squared depth on a normal nearest
edge is then strictly larger than \(f(t)=(4-t^{-4})/t\).
For \(1\le t\le17/16\), \(f(t)\ge3\), which the real covering
radius directly forbids. The function \(f\) has a single interior
maximum (its derivative has numerator \(5-4t^4\)); for
\(17/16\le t\le3/2\), its endpoint values show that \(f(t)>5/2\).
For a normal opposite pair of raw squared difference 4 the
depth must exceed \(g(t)=(4-2t^{-4})/t>9/4\) on that interval;
again the derivative has one maximum and the assertions are
exact at each endpoint.

Now restrict shifts to the half-coordinate torus
\((\tfrac12\mathbb Z^{12})/H\). A phase is encoded by a twelve-bit
mask \(m\) giving its fractional part \(m/2\), followed by a
five-bit integer syndrome \(s\). The encoding is \(32m+s\).
The exact squared depth in quarter-units is
\[
  4q(m,s)=\mathrm{wt}(m)+4
     \min_{\text{syndrome}(A)=s}\#(A\setminus m).
\]
A half coordinate has two minimizing integer translates with
the same cost, so flipping its bit costs zero; all other
unit flips cost one. This reduces to a 32-state subset update.
The group law for adding two phases has the syndrome carry
\(\operatorname{syn}(m\wedge m')\). To subtract a phase
subtract also its half mask; in the binary representation
this is its addition inverse obtained by toggling the
syndrome by \(\operatorname{syn}(m')\).

Every normal FCC octahedron comprises the six points
\(o\pm e_i\) about an odd-sum \(o\), in raw coordinates.
All nonopposite pairs are normal nearest edges. Thus the six
shifts on it must have depth at least \(11/4\) on all edges
and \(5/2\) across the opposite pairs when \(t\le3/2\);
round up the necessary strict thresholds using the quarter
grid. The exact finite replay proves no such six labels
exist. Translate a first vertex to 0. Its opposite can use
4956 labels of depth at least \(5/2\); every other vertex is
a common \(11/4\)-neighbor of that pair. For each such pair
the checker branches on the second opposite pair in this
common neighborhood, and tests the final pair among their
common neighbors. It covers all 827136 middle opposite
pairs, with no allowed final pair in any branch.

Consequently no half-grid assignment on this intact FCC
normal stock can improve the benchmark at \(1\le t\le3/2\).
For an additive assignment (or its translate) this in fact
excludes all \(t\): phases in the half-coordinate torus have
order dividing four, so a fourfold multiple of any nearest
normal edge has zero horizontal shift and squared normal
separation \(16\gamma\). Hence \(\gamma\ge1/4\) and a strict
gain would require \(t<\sqrt2<3/2\).

This leaves larger \(t\) with nonlinear half-grid assignments,
off-grid shifts, deformed normal geometry and incomplete
horizontal layers. As one concrete metric target,
\((t,\gamma)=(8/5,1/8)\) would give about 35% gain if an
assignment existed. A literal 12-phase \(5/2\)-separated
motif is included in the replay; it is only enough to clear
a local neighborhood and is not a global FCC labeling.
At that metric, pairs at raw FCC squared separations
2, 4 and 6 need half-grid depths at least \(5/2\); pairs at
larger raw separation can use shallower differences but must
still obey the exact distance budget until raw separation 64.
The easy additive periodic continuation has already been
excluded above.

Replay: node verify.js. It evaluates all distances and the
combinatorial cases with integer arithmetic.
