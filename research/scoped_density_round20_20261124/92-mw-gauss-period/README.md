# An exact 23-step similarity in the MW92 cusp filtration

Let \(p=47\), \(T\) the order-47 translation of the pinned [Shioda
MW92 parent](../../scoped_density_round17_20261124/92-MW-cusp-jets),
and \(D=1-T\). For **every** translation-invariant full sublattice
\(M\), not just \(L\), \(D^{23}M\) is similar to \(M\) with squared
scale \(47\). The similarity is a single integer operator:
\[
 B=\sum_{j=1}^{46}\left(\frac{j}{47}\right)T^j.
\]
In particular the profile \((a+23,b+23)\) at index \(47^{a+b+46}\)
is in density-preserving bijection with the profile \((a,b)\) at index
\(47^{a+b}\).

Set \(\chi(0)=0\). The coefficients of \(B(T)B(T^{-1})\) are
\(\sum_j\chi(j)\chi(j+k)\): 46 for \(k=0\), and \(-1\) otherwise.
For nonzero \(k\), the terms \(j\ne-k\) equal \(\chi(j/(j+k))\);
this argument of \(\chi\) runs through all residues except 1.
Because \(\sum_jT^j=0\), this product is \(47I\).
Also \(B(T^{-1})=-B(T)\) because \(\chi(-1)=-1\), hence
\(B^2=-47I\). Since \(T\) is a Gram isometry, \(B\) multiplies
every squared norm by 47.

Over the completed discrete valuation ring \(\mathbb Z_{47}[\zeta_{47}]\),
\(D\) is a uniformizer, and \(47\) has valuation 46. Thus \(B\)
has valuation 23 and \(BM=D^{23}M\) at 47. At every other prime
both \(B\) and \(D\) act as units on the lattice, with respective
inverses \(-B/47\) and \((-\sum_j jT^j)/47\). This proves equality
of the integral lattices. Correspondingly \(D^{-23}\) takes every
sublattice with \(a,b\ge23\) back to an integral sublattice of the
parent: it does so at 47 by the elementary-divisor classification,
and at every other prime because \(D\) is a unit. This is the inverse
density-preserving bijection.

In particular \(m(D^{d+23}L)=47m(D^dL)\), for \(m\) the squared
minimum. The period 23 is the smallest positive shift that can give
a similarity of two stem lattices: their covolumes have ratio
\(47^{2e}\) at shift \(e\), so the required squared scale is
\(47^{e/23}\). It must be rational, being a ratio of positive squared
minima of integral lattices. Therefore \(23\mid e\).
This is a symmetry reduction, not a bound on all primitive profiles
or a new density record.

## Literal additional witnesses

The neighboring witness file verifies \(m(D^5L)\le30\) and
\(m(D^8L)\le32\). The latter is **exact**, by the pinned
geometric lower bound \(m(D^dL)\ge16+2d\) at these depths.
The rows came from Gram LLL in rank-46 reflection eigensublattices;
norms and membership are independently checked by integer arithmetic,
without relying on LLL. The norm-32 row was discovered at depth six
and independently passed two more divisibility tests.

Since \(30^{46}<16^{46}47^8\), the first witness applies to all
profiles with \(b\le5\) at index at least \(47^8\); since
\(32^{46}<16^{46}47^9\), the second applies when \(b\le8\) at
index at least \(47^9\). In particular the latter excludes every
mixed elementary-divisor type \((1,8),(2,7),(3,6),(4,5)\) at
index \(47^9\), a total of **517,636,666,560** sublattices, since the
four counts add to \(48(47^6+47^4+47^2+1)\),
without searching those profiles. Any primitive type at that index
is outside this argument.

Run OPENBLAS_NUM_THREADS=1 python verify_gauss.py, optionally
passing the pinned round-17 input.json.gz as the first argument.
The round-17 geometric and height audit provides the lattice
provenance.

## General form of the similarity

For any odd prime \(p\), if an integral Euclidean lattice \(L\) has a
fixed-point-free integral isometry \(T\) of order \(p\), the same
construction works with \(B=\sum_{j=1}^{p-1}\chi(j)T^j\) and shift
\(e_0=(p-1)/2\). The two Gauss-sum identities are
\(B(T)B(T^{-1})=pI\) and \(B(T^{-1})=\chi(-1)B(T)\). Thus
\(B^2=\chi(-1)pI\), and squared norms are multiplied by \(p\),
regardless of that sign. Fixed-point-free and prime-order imply
\(\Phi_p(T)=0\), as needed both for the convolution identity and for
\(D\)-adic valuation. The valuation of \(B\) is \(e_0\); away from
\(p\), both \(B\) and \(D\) remain units. Therefore \(BM=D^{e_0}M\)
for every \(T\)-invariant full sublattice \(M\).

If the rank is \(n\), the characteristic polynomial of \(T\) is
\(\Phi_p^{n/(p-1)}\), so shifting a stem by \(e\) multiplies its
covolume by \(p^{en/(p-1)}\). Any similarity between those two stems
would have squared scale \(p^{2e/(p-1)}\), necessarily rational
because the two squared minima are positive integers. Thus the
half-ramification shift \(e_0\) is minimal.

For an independent geometric verification, polynomial_witnesses.json
also gives two literal sections in (mathbb F_{47}[i]), (i^2=-1).
The coefficients are encoded (a+47b) for (a+bi), in increasing
degree. Run python verify_polynomial.py to check
(Y^2-X^3=t^{47}-t) and degrees (22,33) at depth three and
(32,48) at depth eight. A polynomial section has no finite
zero-section intersection. At infinity its formal parameter
(-z^{-8}X(1/z)/Y(1/z)) has order respectively three or eight, so
the total heights are 22 and 32. The exact finite-fibre
interpolation used to discover the coefficients is not invoked
in this replay; the identities are checked coefficientwise.
