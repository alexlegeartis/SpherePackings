# MW92: complete obstruction at invariant index 47^5

All **234,330,768** translation-invariant index-47^5 sublattices of
[Shioda's MW92](https://doi.org/10.2307/2374791) contain a nonzero
vector of squared norm at most 24. Their center density is therefore
at most $(24/16)^{46}/47^5<0.550$ times the parent. This is a
complete obstruction for this index and invariant family, not for
arbitrary sublattices or arbitrary 92-dimensional packings.

The [previous exact screen](../../scoped_density_round20_20261124/92-mw-cusp-index4-5)
covered all mixed profiles and left 53,016 primitive profiles, over
exactly two three-jet prefix orbits. Here two new literal norm-24
vectors in $K_2=D^2L$ normalize to $(2,27,21)$ and $(6,24,0)$.
Each gives 12 orbit prefixes under the previously checked normalizer,
with 47^2 lifts per prefix. These close exactly the two residual
orbits. Each acceptance in the replay is an integer norm, jet, or
finite orbit verification.

## How the two short sections were found

Let $E$ be one eigenspace of a lattice involution. Hold the
orthogonal component of a known short depth-two row $v$ fixed,
and vary the component in $E$. In the six coordinates of $K_2/K_5$,
membership in a length-three prefix $(a,b,c)$ gives the three
linear congruences
$$
s(w)=(w_1-aw_0, w_3-aw_2-bw_0, w_5-aw_4-bw_2-cw_0)=0
       pmod {47}.
$$
On the chosen rank-46 integral eigensublattice of $K_2$ these
conditions have rank three. The kernel has index 47^3 in the
eigensublattice, so the corresponding Gram determinant is
multiplied by 47^6. A solution $z_0$ to $s(z_0)=-s(v)$ reduces
the search to finding close vectors in the affine lattice
$v+z_0+ker s$. The squared norm splits across the held
orthogonal coordinate and the rank-46 translated form. Reduced
Gram bases and bounded sphere enumeration proposed the two rows,
in the positive coefficient-Frobenius and positive anti-linear
reflection charts, respectively. No completeness of this search
or of any short shell is used by the certificate.

## Beyond this index

The new proof uses only witnesses of depth at least two.
[Non-invariant and higher-index consequences](NONINVARIANT.md)
therefore exclude all 4,985,761 mixed invariant sublattices at
index 47^6, two mixed types at index 47^7, and large explicitly
counted sets of sublattices containing $K_5$ without assuming
translation invariance.

Run `OPENBLAS_NUM_THREADS=1 python -B verify_fifth_complete.py`
using the preceding published proof packets. For a standalone local
copy, give the Round 20 fifth-jet script as the first argument,
followed by the arguments listed in its own instructions.
