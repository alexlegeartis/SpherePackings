# 92D: every translation-invariant index-\(47^3\) sublattice has a short vector

Let \(L\) be Shioda's rank-92 lattice in the cusp coordinates of
[Round 17](../../scoped_density_round17_20261124/92-MW-cusp-jets).
Set \(p=47\), let \(T\) act by translation, and write \(D=I-T\).
All **106,080** translation-invariant index-\(47^3\) sublattices
contain a nonzero vector of squared norm at most 20. Since
\(20^{46}<16^{46}47^3\), none improves the parent center density.
This does not exclude arbitrary non-invariant sublattices or
assert the best minimum attained at this index.

## The quotient problem

In the completed discrete valuation ring
\(\mathbb Z_{47}[\zeta_{47}]\), \(D\) is a uniformizer,
the residue field is \(\mathbb F_{47}\), and the completed
\(L\) is free of rank two. Indeed \(\Phi_{47}(1-X)\) is Eisenstein.
See the local argument of
[Round 18](../../scoped_density_round18_20261124/92-mw-cusp-index2).
Every translation-invariant quotient of order \(p^3\) is
annihilated by \(D^3\).
There are 48 profiles of elementary-divisor type \((1,2)\),
one for each line between the first and second jet. They all
contain \(D^2L\) and hence its checked norm-20 vector.
The remaining type \((0,3)\) profiles are primitive projective
lines over \(\mathbb F_p[D]/D^3\). Their count is
\((p+1)p^2=106032\): normalize the first unit coordinate to 1,
then the other is an arbitrary three-digit element if the first
coordinate is a unit, or a two-digit multiple of \(D\) otherwise.

Build the jet matrix \(J_3\) on basis rows using \(A,H,R\)
from Round 17: emit \(xA\bmod p\), subtract the lifted residue
section via \(R\), multiply the result by \(H/p\) for \(H=pD^{-1}\),
then repeat. \(J_3\) has rank six, and \(D\) shifts its successive
coordinate pairs. Its kernel is \(D^3L\) by the index.
The integral isometries \(U,C,F\) normalize \(T\); the replay
verifies their induced jet maps. A short vector of depth \(d\le2\)
gives a normalized ratio modulo \(D^{3-d}\), witnessing all
\(p^d\) primitive lines with that prefix.

## Short-vector coverage

In the chart \(Y^2=X^3+t^{47}-t\), inversion \(t\mapsto-1/t\)
acts on sections by \(X'(t)=t^{16}X(-1/t)\), \(Y'(t)=t^{24}Y(-1/t)\).
The replay checks the equation on the 92 pinned polynomial
sections and reconstructs its integer matrix \(S\) by their
exact height pairings. Solving modulo 65521 only proposes
entries: the full pairing identity, \(SGS^\mathsf{T}=G\),
and \(S^2=I\) are checked over the integers.

Start from 119 frozen rows: 17 Round 17 norm-16 rows, two
norm-18 rows, six Round 18 witnesses, 92 independent norm-16
rows, and two depth-two norm-20 rows. Apply \(T^jS\) to each,
\(0\le j<47\). Then take both signs of the sum of one frozen
row and another transformed by the identity or
\(MT^jS\), \(M\in\{I,U,U^2,U^3,C,C^2,F\}\).
Use the actual Gram bilinear form to retain only norms 16,
18 and 20; after deduplication up to sign there are 1,460
verified short rows. Apply \(T^jS\) to these too and close
their normalized prefixes under \(U,C,F\).

The depth-two, depth-one, and depth-zero prefix counts are
12, 558, and 106,002, giving 26,508, 12,972, and 66,528
incremental exclusions after overlap, respectively.
Adding the 48 mixed profiles gives 106,056. One literal
norm-18 row from the earlier cusp search hits 12 more through
the verified normalizer. The final 12 lines form one orbit,
with representative \((1,3+12D+16D^2)\). One further
literal norm-20 row has prefix \((1,3+12D)\) at depth one
and therefore witnesses that representative. Its verified
normalizer images witness the whole remaining orbit.
Both literals are in the checker; the discovery archives
are not needed.

Run OPENBLAS_NUM_THREADS=1 python -B verify.py (NumPy).
It uses the preceding exact Round 17 and 18 packets;
their paths may be supplied explicitly. All accepting
decisions use exact lattice vectors and quotient arithmetic.
The [index-\(47^4\) extension](FOURTH.md) covers part of
the next index and leaves the others open.
