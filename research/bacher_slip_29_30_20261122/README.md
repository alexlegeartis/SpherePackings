# Shifted half-phase lattice packings in 29D and 30D

These are **lattice-only** advances. The stronger nonlattice B29* and T30
bounds in [Cohn’s packing table](https://cohn.mit.edu/sphere-packing/)
remain unchanged.

| Dimension | Previous catalogue lattice center density | New center density | Relative gain |
|---:|---:|---:|---:|
| 29 | 1/√3 (B29) | **4/√33** | **20.6045%** |
| 30 | 3^(27/2)/2^22 (Q30) | **8/√105** | **18.5819%** |

The named comparisons are from the [Nebe–Sloane catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/density.html).
The 28D parent is Bacher's bac28aR in [Martinet's exact Gram collection](https://jamartin.perso.math.cnrs.fr/Lattices/bacher.gp),
from [Dense lattices in dimensions 27–29](https://doi.org/10.1007/s002220050180).
This comparison does not certify worldwide priority or assert optimality.

## Construction

Take the 28 horizontal basis generators of the given Gram G, with
zero vertical component. From the binary masks in witness.json,
form their 0–1 coefficient vectors a,b in little-endian order.
Append one generator with horizontal coefficient vector a/2+d/12
in 29D. Append two generators with horizontal coefficients
a/2+d/12 and b/2+d/12 in 30D. The listed integer d has squared
norm 24. The vertical Gram is [11/6] in 29D and
[[11/6, -2/3],[-2/3,11/6]] in 30D.

The shifts d/12 are essential. The squared clearance of each original
halved phase is 6, but its nearest vectors span a hyperplane. Moving
off that hyperplane increases the clearance, permitting a lower
height. Checking the shift also checks the even vertical layers,
which are no longer horizontally unshifted.

The horizontal lattice has squared minimum 8. For a nonzero
vertical coefficient vector, if a coefficient has absolute value
at least 3, the smallest eigenvalue of its vertical Gram alone
proves squared norm greater than 8. In the remaining sign-symmetric coordinate
range, the verifier enumerates every vertical vector below 8 and
checks a finite exact shifted CVP exclusion at the appropriate
horizontal threshold. This accounts for all lattice translates.
Both extended lattices therefore have and attain squared minimum 8.

Schur complementation gives Gram determinants 1,107,296,256 and
1,761,607,680. With squared radius 2, the exact center-density
squares are 16/33 and 64/105.

## Exact replay

Run python3 verify.py using Python 3.9+ and a C++17 compiler.
The Python code verifies the catalogue Gram hash, derives a rational
LDL factorization, and generates all shifted CVP inputs. In the DFS,
each LDL term is rounded down at integer scale 2^20. The branch bound
uses the exact integer inequality for that rounded term, so it
cannot exclude a truly short vector. Integer square-root estimates
are corrected in both directions and every leaf is checked against
the literal input Gram. A separate independent audit uses a different
unimodular basis and integer LDL enumeration.
