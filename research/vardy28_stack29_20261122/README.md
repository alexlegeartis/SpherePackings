# A 29-dimensional nonlattice packing of center density 1

A periodic packing with exact center density **1** and occupied density
`4.82878227389174e-5…`. This is a **41.421356%** improvement over the
29-dimensional benchmark `1/sqrt(2)` in [Cohn's unrestricted table](https://cohn.mit.edu/sphere-packing/#data),
the 29-dimensional packing of [Vardy (1999)](https://doi.org/10.1090/S0002-9947-99-02169-8).
See also the [Nebe–Sloane catalog](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/density.html).
It is not a lattice construction or a claim of optimality.

## Construction

Use two binary codes of length 28. Split their coordinates into two halves of
14. Write (A) for the linear map on 14-bit masks whose images of successive
unit vectors are the successive cyclic left rotations of the 14-bit mask
**471**. Put (C={(x,Ax)}) and (D={(A^Ty,y)}). Each is even, has
rank 14 and minimum Hamming weight 8, and the two are dual under the binary
dot product. These facts are checked directly. Their density-doubled stock
uses the form introduced by Vardy in dimension 28.

Let (D_{28}={zinmathbb Z^{28}:sum z_iequiv0pmod2}), independently
of the binary code (D). The horizontal point set is
[
P=igcup_{cin C}(-2c+4D_{28}) cup  igcup_{din D}(mathbf 1-2d+4e_1+4D_{28}).
]
Choose (b=(715,15538)) in the same two-half bit notation. Within each
half the least significant bit is the first coordinate. Its minimum Hamming
distance is exactly 6 to **each** of (C) and (D). Put an unshifted copy
of (P) at height (2jsqrt 8) and (P+2b) at height
((2j+1)sqrt8), for every integer (j).

## Separation and covolume

* Within the same layer and parity stock, distinct binary cosets force at
  least eight coordinates of absolute difference at least 2. A nonzero
  pure period in (4D_{28}) has squared length at least 32.
* Opposite-stock pairs within a layer have all coordinates odd. To obtain
  absolute values all 1 would require carries of parity
  (operatorname{wt}(c)-ccdot d=0pmod2): a carry 1 is needed exactly
  when (c_i=1,d_i=0). The (4e_1) offset requires odd carry parity.
  Thus at least one coordinate has absolute value at least 3, for
  squared distance at least 36.
* In adjacent layers, same-stock horizontal differences have at least six
  coordinates equal to 2 modulo 4, by the simultaneous distance test and
  linearity. Their horizontal squared distance is at least 24. Opposite
  parity forces horizontal squared distance at least 28. Add 8 from the
  vertical separation in either case. Different layers farther apart
  already have vertical squared distance at least 32.
* The attained squared minimum is 32, for instance along the period
  (4(e_1+e_2)) in an unshifted layer. The horizontal period determinant
  is (2cdot4^{28}), with exactly (2^{15}) horizontal points per
  period, distinct modulo 4. The vertical period is (2sqrt8), doubling
  the multiplicity. The volume per point is (2^{42}sqrt8). The
  radius is (sqrt8), so the center density is
  ((sqrt8)^{29}/(2^{42}sqrt8)=1).

The set contains zero but does not close under addition. In the
height-zero layer, take the even point (-2(1,471)) and the odd point
(mathbf 1+4e_1). Their sum would require ((1,471)in D) by its
mod-4 residue. The verifier confirms that this word is in (Csetminus D).

Run `python verify.py`; only the Python standard library is needed.
The verification enumerates all words of both codes from the 14-bit
generator, checks the two depths and parity properties, and checks the
exact normalization. The separation argument covers every pair of
periodic translates, not just those in one fundamental domain.
