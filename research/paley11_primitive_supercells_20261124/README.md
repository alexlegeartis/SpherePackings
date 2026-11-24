# Two 11-dimensional packings tying the known center density with larger primitive cells

The files give explicit motif points in coordinates doubled for storage. Use
period lattices `2D4 × 2Z^7` for `motif144.txt`, and
`2D4 × 2D4 × 2Z^3` for `motif288.txt`. Here
`D4 = {x in Z^4 : sum(x) is even}`, with the usual Euclidean metric. The
construction pools block-triality images of lifts of the 72-site Paley motif
and selects a compatible motif from the pool. To obtain the three orientations
on a block, use the identity, the normalized 4-by-4 Hadamard transform, and
a one-coordinate sign change composed with that transform. These preserve
D4: the Hadamard transform is orthogonal, outputs are integral on an
even-sum integer input, and their sum is twice the first input coordinate. The literal files are sufficient
to reproduce the packings; no optimizer is required.

| Sites | Minimum squared distance | Primitive period covolume | Center density |
|---:|---:|---:|---:|
| 144 | 4 | 4096 | 9/256 |
| 288 | 4 | 8192 | 9/256 |

Run `python3 verify.py` (standard library only). It checks every pair by exact
nearest-period reduction; each four-coordinate block is projected to the
even-sum integer lattice after dividing the doubled displacement by four. A
nearest parity-correct lift differs from the coordinatewise nearest lift in
at most one coordinate. The script also checks every possible extra period:
an extra translation modulo the listed lattice must take the first site to
another motif site. Only the identity shift works.

Primitive period covolume is invariant under Euclidean congruence. Therefore
the two examples differ from each other and from any packing of this same
minimum having period `2Z^11` (whose primitive covolume is at most 2048).
They give constructive starting points with more site coordinates for a
search without imposing the smaller period, not a density increase or a
claim of historical priority. The density
comparator is Conway and Sloane's
[11-dimensional construction](https://neilsloane.com/doc/Me185.pdf).
