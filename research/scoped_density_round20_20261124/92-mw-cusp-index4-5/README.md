# MW92: every translation-invariant index-47⁴ sublattice falls short

Let \(L\) be [Shioda's rank-92 Mordell–Weil lattice](https://doi.org/10.2307/2374791),
of squared minimum 16 and determinant \(47^{14}\). In the cusp chart
\(Y^2=X^3+t^{47}-t\), let \(T\) be translation by one and \(D=I-T\).
**All 4,985,761 translation-invariant index-47⁴ sublattices have squared
minimum at most 22**. Each has center density at most
\((22/16)^{46}/47^4<0.472\) times the parent. This is a finite family
upper bound, not a bound on arbitrary sublattices or general packings.
No attainment at minimum 22 is claimed.

There is also a literal norm-22 vector in \(D^3L\). The geometric lower
bound from [Round 17](../../scoped_density_round17_20261124/92-MW-cusp-jets)
therefore gives the **exact squared minimum 22** at depth three.
The literal came from Gram reduction of a rank-46 reflection
eigensublattice. The replay uses the vector itself, not reduction output.

## Exhaustive fourth-jet witness count

The completed translation ring at 47 is the DVR
\(\mathbb Z_{47}[\zeta_{47}]\) with uniformizer \(D\), and \(L\)
is free of rank two over it. Modulo \(D^k\) for \(k\le5\), the
ring is \(\mathbb F_{47}[D]/D^k\). The integer jet compiler
verifies the quotient dimension \(2k\) and that \(D\) shifts
successive coordinate pairs. This identifies the first \(k\)
digits with any choice of two lifts of a basis of \(L/DL\).
The finite isometries \(U,C,F\) normalize \(T\); the checker
verifies their induced maps on all 92 integral coordinate basis rows.

There are three index-47⁴ types. The 2,256 type-(1,3) and the
one type-(2,2) sublattices all contain \(D^3L\). Type (0,4)
gives 4,983,504 projective lines over the four-jet ring.
To normalize, set the first unit coordinate to 1, or use the second
coordinate if the first is a nonunit. A vector of depth \(d\)
then determines a ratio modulo \(D^{4-d}\) and belongs to every
one of the \(47^d\) lifts of that ratio. Keep shortest witnessed
prefixes first, so the contributions are disjoint.

The [Round 19 screen](../../scoped_density_round19_20261124/92-mw-cusp-index3/FOURTH.md)
had reached 3,665,047 of these sublattices. Appending \(MT^j\)
to its \(MT^jS\) pair operators, for the same seven \(M\) and 47 shifts,
yields 11,876 pair witnesses up to sign from the original rows alone.
The new depth-three row and its normalizer orbit then cover six full first
digits and every mixed profile. The 27 further literal rows, in three
discovery batches of 9, 6, and 12, are directly checked over the
integers. The first nine are depth two, from pairing norm-18
vectors at a common cusp and matching two-jet residues.
The last twelve cover the last twelve primitive orbits, via independent
translation of the two input sections before pairing. Only their
resulting witness rows, not their discovery archives, are needed.

The primitive counts from shortest-prefix depths 3, 2, 1, 0 are
622,938; 4,334,058; 12,972; and 13,536. Their sum is 4,983,504.
With the mixed types, this exhausts all 4,985,761.

## Fifth jet: a finite partial screen

At index \(47^5\), \(24^{46}<16^{46}47^5<26^{46}\), so all
norm-24 witnesses also obstruct density improvement. Every mixed
profile of type (1,4) or (2,3) contains the pinned literal norm-24
row in \(D^4L\); this covers all 106,080 mixed profiles.
In the primitive type, 12 witnessed leading directions by themselves
cover \(12\cdot47^4=58,556,172\) profiles; the new orbit at
slope 46 is disjoint from the original six. Direct depth-three and depth-two rows give 57,310,296 and 118,305,204
additional primitive cases after overlap, respectively. Together with
the mixed cases, the exact direct screen covers **234,277,752 of
234,330,768 profiles**. This uses literal rows from common-cusp signed
pairs and three-row combinations with matching second-jet residues,
not a full shell enumeration. The remaining 53,016 primitive cases lie
over only 24 three-digit prefixes, in two sets of 12 related by the
verified normalizer. Representatives, written as five pairs of jet
coordinates, are `(0,0; 0,0; 1,2; 0,27; 0,21)` and
`(0,0; 0,0; 1,6; 0,24; 0,0)`. No statement follows about these
remaining profiles.

Run from this folder with NumPy:
```sh
OPENBLAS_NUM_THREADS=1 python -B verify_fourth.py
OPENBLAS_NUM_THREADS=1 python -B verify_fifth.py
```
Optional explicit arguments are Round 19's index-three verify.py,
Round 17's cusp-jet verify.py, and Round 18's index-two directory,
in that order. Both scripts use the exact earlier proof and freeze the
Gram, isometries, and translation to that input. The fifth extension uses only direct prefixes through length three;
it does not enumerate all quarter-billion profiles individually.
