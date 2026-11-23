# A six-translate 59-dimensional sphere packing from a matched Ne64 normal

The packing fraction is **2.64616393271302704634577 × 10⁻¹²**,
an improvement of **0.32458801879956541093%** over the repository's
[previous six-translate 59D packing](../../constructions/dim59/).
The center density is **275413.65717634773482453485**. This is an
explicit periodic nonlattice packing, not a lattice-only improvement
or an upper bound.

The two constructions differ in the shape of the normal lattice. The
previous one used an off-diagonal path of four −3 entries in its normal
Gram numerator. Here the special edges form a matching: there are two
−4 entries on disjoint edges, and all remaining off-diagonal entries
are −2. The new projected cluster has a slightly greater squared
diameter, but the smaller period covolume more than compensates for it.

## Exact construction and inherited inputs

The only inherited ambient facts are **minimum squared length 12** and
**3-modularity** of the [Nebe–Sloane catalogue lattice Ne64 (CQ64b)](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ64b.html).
See the MINIMAL_NORM and PROPERTIES entries and Nebe's references there.
With catalogue Gram G and integral unimodular change T from change.txt,
put C=T G Tᵀ. The alternative ambient Gram is J=3 C⁻¹, another coordinate
model of Ne64 by 3-modularity. The literal J is in data/ambient_gram.txt.

The rows of A=normal.txt are coefficients in the **dual basis** of J,
so their normal Gram is K=A J⁻¹ Aᵀ=A C Aᵀ/3. Its exact numerator is

    12  -2  -2  -2  -2
    -2  12  -4  -2  -2
    -2  -4  12  -2  -2
    -2  -2  -2  12  -4
    -2  -2  -2  -4  12

and det K=32768/81. The minor on columns 0,1,3,5,25 is unimodular,
which proves that the normal is primitive and every integer normal
label is attained by an ambient point. Choose the six labels
0,e₁,e₂,e₃,e₄,e₅. Their squared normal distances are 3/8 or 7/16.

The period basis N is the entire integer kernel of A, as verified by
its identity block on the complementary columns and the unit minor.
Its Gram Q=N J Nᵀ has determinant 32768·3²⁸. The six phase shifts,
with denominator 96 in that basis, are obtained by projecting the
corresponding ambient lifts perpendicular to the five normal rows.
Explicit period Gram and shift numerators are in data/.

A difference between centers lifts to a nonzero ambient vector with
normal component of squared norm at most 7/16. Its horizontal squared
norm is at least 12−7/16=185/16. The projected translates are distinct:
if two coincided, their lift difference modulo the period would be a
nonzero purely normal ambient vector of squared norm below 12.
The ambient norm-12 contact in contact.txt has normal label −e₄ and
attains the diameter 7/16, proving exact minimum 185/16. The
construction is not a lattice: twice the first nonzero motif shift
belongs to none of the six classes modulo the period, checked exactly.

Hence the squared center density is
`36 (185/64)^59 / (32768·3^28)`. Both checkers compare it as an exact
rational to the predecessor's `3^33 401^59 / (2^296 25 13^61)`.
Ordinary packing fraction is obtained by multiplying the center
density by π^(59/2)/Γ(61/2).

## Replay

  * Primary exact: python verify.py. Requires SymPy. To regenerate the
    literal matrix and six shifts, run python verify.py data.
  * Independent exact: python verify_independent.py. Uses the Python
    standard library and verified literal matrix relations, a separate
    integer determinant routine, exact rational shift distances and the
    literal shortest contact. It neither imports nor calls verify.py.

The archive starts with the literal normal matrix already selected.
No numerical exploration, optimization, or ambient minimum enumeration
is used in either proof. The construction is not a global-priority claim
beyond the compared and specified predecessor.
