# A three-point cap on heavy radical packets in the 12D equal-scale mask

Work in the equal-scale \(E_6\times E_6\) ternary motif mask:
\(V=(E_6/3E_6)\), and distinct motif sites must have sum of the two
least squared lift norms at least 8. Let \(R\) be the order-three radical
of the mod-three bilinear form, and call a coset of \(R\times R\) a
packet. A packet is heavy if it contains at least five of its nine
possible sites.

**There are at most 108 heavy packets.** In particular, a motif whose
every occupied packet is heavy has at most \(108\cdot9=972\) sites,
down from 1215, and so does a union of full packets. This does not
exclude arbitrary motifs with the 923 sites needed to beat
Coxeter–Todd. It is not an unrestricted sphere-packing bound.

## Reduction to the five-dimensional quadratic space

For \(W=V/R\), assign types \(0,1,2,3\) to its zero, the two nonzero
quadratic values, and nonzero isotropic elements, respectively.
The minimum squared lift norms among three lifts are \(0,2,4,6\);
see the exact enumeration in
[Round 18](../../scoped_density_round18_20261124/12-e6-ternary/README.md).
If the two quotient difference types of distinct packets are \(i,j\)
with \(i+j<4\), some fixed radical offset gives a forbidden difference.
The phase sets in two heavy packets must realize every fixed offset:
two subsets of a group of size nine with total size greater than nine
intersect under every relative translation. Hence the occupied heavy
packet labels form a code on \(W\times W\) whose nonzero differences
have \(i+j\ge4\). For distinct packets neither type can then be zero.

We exclude a code of exactly 109 labels, which also excludes every
larger one by taking a subset. This uses actual joint triangles,
not saturation of low-occupancy packets and not just scalar
autocorrelation on \(V\times V\).

## Exact finite three-point certificate

Use residue coordinates \((0,w)\) in \(E_6/3E_6\), which are unique
representatives of \(W\) since the radical generator's first coordinate
is nonzero. For ordered \((x,y)\) in \(W\), its stabilizer orbital has
labels \(\operatorname{typ}(x),\operatorname{typ}(y),
\operatorname{typ}(x-y),\mathbf1_{x=-y}\). These determine the binary
Gram and all linear dependence data (zero and \(x=y\) are already
detected); on a nondegenerate quadratic space the subspace isometry
extends to the whole space by Witt's extension lemma. There are 40
orbitals. The verifier checks the six vertex permutations by direct
coordinate subtraction on all ordered pairs. Take the product of
these orbitals for \(W\times W\) and symmetrize on swapping the factors
and the three triangle vertices. Excluding all nonpacking
differences leaves 56 nonnegative variables, including coincidences.

For a putative code \(C\) of size \(N=109\), the entry \(t(x,y)\)
of the rooted moment matrix is the average under the stabilizer
of the origin of
\[
  \frac1N\sum_{z\in C}1_C(z+x)1_C(z+y).
\]
It depends only on the specified orbitals and is unchanged when the
triangle vertices are permuted, by reindexing the translating root
within \(C\). The analogous entry \(a(x-y)\) with translating root
\(z\) unconstrained is the symmetrized pair autocorrelation, and
equals \(t(0,x-y)\). Thus both the matrix \(t(x,y)\) and its complement
\(a(x-y)-t(x,y)\) are positive semidefinite: they are nonnegative
averages of outer products of translated code indicators, over roots
inside or outside \(C\), respectively. Also
\(t(0,0)=1\), and for each fixed \(x\),
\(\sum_y t(x,y)=N\,t(0,x)\).
The verifier sums the latter identities by the 16 product shell types
using exact orbital multiplicities.

The supplied certificate uses 22 positive rational multiples of
Rayleigh inequalities from these two matrices and unrestricted
rational multiples of the row-sum equalities and \(t(0,0)=1\).
It has right-hand constant 1, while every coefficient of a
nonnegative triangle variable is strictly negative. This is the
contradiction. The Rayleigh vectors are encoded as integer linear
combinations of products of ordinary 243-entry integer vectors;
no spectral approximation or SDP solver is used by the verifier.

Replay: python -B verify_packets.py (NumPy for exact bounded integer
arrays, and Python integers and rationals for the separating check).

## Independent local observation: at most three sites in an actual fibre

Fixing the entire first coordinate in \(V\), the second coordinates in
that fibre have all differences of depth at least 8. Translate one to
zero; its remaining points lie among the 146 depth-\(\ge8\) residues.
The graph of allowed pairs on these residues is triangle-free, so each
actual coordinate fibre contains at most three sites. The exact
replay also counts the anchored edges: 1440 noncollinear triples of
depths \((8,8,8)\), 216 noncollinear triples of depths \((8,8,12)\),
72 collinear triples of depths \((8,8,8)\), and one radical line of
depths \((12,12,12)\). This is separate from the packet bound.
Replay: python -B verify_fibres.py.

## Separation from ordinary two-point tests

On the \(4\times4\) quotient shells, per-displacement weights
\(p_{00}=1,\ p_{13}=p_{31}=1/160,\ p_{22}=1/150,\ p_{33}=1/800\)
and zero elsewhere have nonnegative Fourier transform and formal
size 135. The checker gives the exact four quotient character rows
and verifies positivity and the size. Filling all nine phase positions
at each formal packet yields the usual two-point pseudo-size 1215
while respecting valid packet local shapes. The three-point proof
identifies a genuine obstruction missing from this two-point model,
not a discrepancy due to tolerance or to packet saturation.
