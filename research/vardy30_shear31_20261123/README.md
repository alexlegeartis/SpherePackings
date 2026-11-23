# An order-four shear packing in 31 dimensions

This explicit periodic packing has squared minimum **32**, exact center
density **sqrt(2)**, and occupied density at least **1.38411060495194e-5**.
This is **16.9233029%** above the listed Q31 lower bound, whose exact center
density is 3^15/2^(47/2). See the [Nebe–Sloane catalogue](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html) and
[Cohn's sphere-packing tables](https://cohn.mit.edu/sphere-packing/#data).
It is not a lattice and does not change the lattice-only lower bound.

## Construction

Split the horizontal coordinates into two blocks of 15. Let A have as its
columns the successive cyclic left rotations of 1491 in 15 bits. Use
the binary codes C={(x,Ax)} and D={(A^T y,y)}. They are even, dual, and
have minimum Hamming weight eight. The binary word b=(32767,0) has odd
weight and Hamming distance seven to both codes. Both assertions are
verified by enumerating 2^15 words in each code.

Let D30 denote the integer vectors of even coordinate sum and e1 the
first coordinate vector. In horizontal space take the [Vardy double-stock
packing](https://doi.org/10.1090/S0002-9947-99-02169-8)
```
P = union over c in C of (-2c+4D30)
  ∪ union over d in D of (1-2d+4e1+4D30).
```
The 31D packing consists of **P + 2jb at height 2j**, for every integer
j. The shear has order four modulo 4D30. In particular the third layer
is not identical to the first: its horizontal shift 4b is the nonzero
carry class, which is essential for this height spacing.

## Separation of every pair of translates

* **Same layer, same stock.** A nonzero binary difference costs at least
  eight coordinates of squared cost four. If it vanishes, a nonzero
  period vector has squared norm at least 32.
* **Same layer, opposite stocks.** All horizontal coordinates are odd.
  For every coordinate to be ±1 the sum of the required carries, modulo
  two, is wt(c)-c·d=0 by evenness and duality. The prescribed offset
  4e1 instead imposes odd carry parity. So at least one coordinate has
  absolute value at least three, and the squared distance is at least 38.
* **One layer apart.** Within one stock, the distance-seven translates
  cost squared horizontal distance at least 28. Between stocks the
  difference remains all odd and costs at least 30. The height adds four.
* **Two layers apart.** Within one stock, different binary residues
  cost at least 32. Equal residues differ by 4b modulo 4D30, so cost
  at least 16 by oddness of b. Between stocks the difference costs at
  least 30. The height adds 16.
* **Three or more layers apart.** The height alone gives at least 36.

The vector (4(e1+e2),0) attains squared separation 32. For the period
4D30 × 8Z the covolume is 2^64; there are four height layers, each
containing 2 × 2^15 distinct horizontal classes. Hence volume per center
is 2^46 and the exact center density is 8^(31/2)/2^46 = sqrt(2).
The height-zero slice is not closed under addition: summing the odd
class of d=0 with the even class of c=(1,1491) would require c in D.
A lattice's height-zero slice would be a sublattice, so the packing is
genuinely nonlattice.

## Replay

```sh
python3 verify.py
python3 independent_audit.py
```

The primary verifier uses exact integer code calculations and rational
Machin bounds to certify the displayed occupied-density floor. The
independent replay enumerates both codes by a separate implementation
and checks the carry, period, density, and nonlattice distinctions.
