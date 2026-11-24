# Exact low-norm obstructions for the quartic p43 84D glues

**Research notes, not a new sphere-packing bound.** This packet rules out
the hoped-for squared minimum 12 on all 32 order-21-invariant Gaussian
Borel glues of a specified quartic Mordell-Weil lattice. It also proves
two scope-specific obstructions to recovering minimum 12 by taking
cosets of finite quotients of the canonical glue. A single standard-library
Python replay verifies the frozen matrices, literal vector memberships,
and finite-field linear algebra.

Run

    python -B replay.py

The files receipt.00, receipt.01, ... are ordered ASCII pieces of a
deterministic gzip-compressed, base85-encoded JSON receipt. They include
the literal integral Gram matrices, isometries, bases, and vectors.
The verifier starts from the previously reconstructed quartic-height
Gram and does not independently rebuild that Gram from elliptic
sections. No solver, network access, binary dependencies, or external
source file is needed for the tests here.

## All 32 Gaussian Borel glues are obstructed

In the pinned twice-height normalization the underlying rank-84
quartic parent has Gram A of determinant $43^{20}2^{42}$, with a
20-dimensional radical modulo 43. The 43-primary discriminant is thus
elementary; identify it with the radical via $u\mapsto u/43$.
The pinned integral order-21 isometry R and Gaussian complex structure
J act on this space. R has ten distinct eigenvalues, with
two-dimensional eigenspaces over $\mathbf F_{43}$. Each eigenspace
is J-invariant and has no proper nonzero J-invariant subspace, since
$J^2=-I$ and $-1$ is nonsquare mod 43.

The discriminant pairing is nondegenerate between inverse-eigenvalue
partners; other eigenpairs pair trivially because R preserves the form.
Therefore an R- and J-invariant maximal (dimension-10) isotropic space
chooses precisely one eigenspace from each of the five inverse pairs.
This produces exactly $2^5=32$ candidates in this scope. Their
overlattices consist of $u/43$ with $u$ integral and residue $u\bmod43$
in the specified chosen eigenspaces. The replay checks even integrality
on the generating residues, verifies the chosen masks, and checks a
literal squared norm at most 10 in every one. This is sufficient for
exclusion without needing any assertion about shell completeness.

The previous strongest pilot candidates #8, #13 and #28 each have a
literal squared norm 8. The four survivors of initial reflected-basis
reduction, #6, #19, #20 and #29, each have literal norm 10. An additional
exact enumeration of the positive fixed lattice of F and of FJ for
those four variants completed in fewer than 34 million nodes per half,
giving respectively 462, 210, 420, 420 norm-10 vectors **per half**
and no norm <=8 in those halves. This is not a global lower bound on
their full rank-84 minima. Other 43-isotropic choices which are not
R- and J-invariant are not excluded.

## Every R-invariant additive quotient is obstructed for the canonical glue

The canonical even Gram G, determinant $2^{42}$, is related to the
raw Gram by the literal numerator basis B with $BA B^T=43^2G$.
The receipt checks its residue space is precisely the canonical mask.
In its integral coordinates we have two literal exact squared-norm-10
vectors x, each with a literal integral ancestor y satisfying

    x = y (R-I).

Thus **every R-invariant additive character into any abelian group**
vanishes on a norm-10 vector: $\chi(x)=\chi(yR)-\chi(y)=0$.
A periodic packing retaining full cosets of that kernel cannot escape
this vector on a retained layer. In particular this rules out every
horizontal rank-four R-fixed rational/integer projection and all
mod-33 candidates in that specific invariance scope. The same argument
does not require a prime modulus or coprimality with 21.

For clarity, trace zero alone would only prove the modular statement
for moduli coprime to 21, because summing all 21 conjugates gives
$21\chi(x)=0$. It would still prove the real/rational fixed projection
statement. The literal integral ancestor is the additional argument
removing that modulus restriction.

The replay independently checks actual q10 witnesses on all five
difference labels of the earlier modulus-33 motif, both index-17
kernels, and explicit index-19 and index-23 kernels, including the
index-19 nontrivial R-eigencharacter (checked by a literal witness, not
by the R-invariance theorem). Any prior 85D lift with one of those
full horizontal kernel/phase layers inherits a blocker at height zero.

## A smaller order-7 invariance still rules out rank-four binary quotients

For $S=R^3$ of order 7, let $P=I+S+\cdots+S^6$ over $\mathbf F_2$.
Every S-invariant binary linear map factors through P, because
$\chi(xP)=7\chi(x)=\chi(x)$. The replay supplies 31 literal q10
vectors whose P-projections are exactly the nonzero elements of a
five-dimensional binary subspace. It checks incidences in all 84
coordinates, not just the 12-column labels used for discovery.
Consequently any such map of rank at most four vanishes on at least
one of those q10 vectors. This does **not** rule out rank >=5,
non-invariant binary maps, or maps invariant only under R^7 (order 3).

## Exact short-vector reduction behind discovery

The real isometries F and FJ of the canonical lattice are involutions
with saturated positive and negative eigensublattices of rank42.
J exchanges the two halves and preserves their Gram. The verifier
checks saturation using an explicit unimodular extension, and each
half Gram has determinant $2^{63}$. A full vector x gives
$u=x+xF$, $v=x-xF$, with $q(u)+q(v)=4q(x)$. Recovery requires
$u+v$ even in the original lattice, not just that the two halves have
small norms. There is an identical FJ version.

The canonical order-43 isometry T commutes with J and satisfies
$FTF=T^{-1}$. Over the Gaussian field its primitive order-43 block
is irreducible of dimension 42. In the unitary eigenbasis F acts
anti-linearly within each one-dimensional eigenspace. Accordingly for
$c_k=\langle x,xFT^k\rangle_{\mathbb C}$, Parseval gives
$\sum_{k=0}^{42}|c_k|^2\ge(43/42)q(x)^2$. For q10 this exceeds 86;
the real and imaginary parts are integers, so some coordinate has
absolute value at least 2. Up to conjugating the involution by T
and exchanging signs by J, one of the F/FJ half-components therefore
has squared norm at most 16 (or vanishes, making x genuinely fixed
in the complementary half). The q<=8 case follows already by
averaging. This small-half reduction is not an exhaustive global
enumeration; it was used to isolate cheap pure-half witnesses.
The canonical half enumerations completed at 32,251,628 and
31,234,779 nodes with 210 vectors of norm10 each, no smaller vector
in those halves. The F/FJ norm10 stocks gave ten disjoint
3612-vector isometry orbits beyond the initial one in the parent.

To reproduce the optional canonical half-node counts:

    python -B verify.py --enum-input
    g++ -O3 -std=c++17 enum_half.cpp -o enum_half
    ./enum_half 550000000 < F_10.in
    ./enum_half 550000000 < FJ_10.in

Inputs are generated with exact rational LDL factors rounded downward.
The search subtracts a conservative floor error from each absolute
center discrepancy before applying its pruning bound, so retained
leaves contain every vector within the specified radius; exact integral
norms at leaves decide the result. The nonzero sign symmetry doubles
counts. The default replay verifies all negative packing implications
directly from literal witnesses, independently of these completeness
counts or discovery tools.
