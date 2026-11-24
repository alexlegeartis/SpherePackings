# Product-sphere bound on the maximum-dot layer of the 96D golden shell

Let L be an even unimodular lattice in dimension 96, with minimum
at least 8 and a self-adjoint operator S with square 5 such that
U=(1+S)/2 preserves L. Assume q(v)=8 implies s(v)=0, where
q(v)=v.v and s(v)=v.Sv. Put X={x in L: q(x)=12, s(x)=0}.
Then every x in X has **fewer than 11,040,000 vectors y in X with
x.y=8**, a factor exceeding 28 below the previous 313 million bound.
The same is true for dot -8 by negation. This is a local theorem,
not an exclusion of the 96D endpoint and not a new packing record.

As in the prior norm-12 local argument, map y injectively to r=x-y.
Then q(r)=8, s(r)=0, x.r=4, and x.Sr=0. The eigenspaces of S
have dimension 48 each. In each eigenspace the projection of r has
squared norm 4, that of x has squared norm 6, and their dot is 2.
Subtract x/3 from r and normalize separately to unit length in the
two eigenspaces, leaving each factor's dimension 47 and pre-normalization
squared norm 10/3. For two distinct r,t, put p=r.t and b=r.St.
Their resulting two inner products are
 (3p-4+3b/sqrt(5))/20 and (3p-4-3b/sqrt(5))/20.
The integer p lies in [-4,4] since r,t are norm 8, distinct, and
not opposites (both dot x=4). We have |b|<=8 sqrt(5)<20 by Cauchy,
and p=b mod 2 since r.Ut is integral. The kernel below is symmetric
in the two factors, so test b>=0 by replacing S with -S, which
preserves all lattice hypotheses.

For each candidate p,b the replay builds the Gram of x,Ux,r,Ur,t,Ut
entirely from U^2=U+1 and self-adjointness. It tests all integral
combinations with coefficients in [-2,2]. Combinations of zero norm
are allowed; every other combination has q>=8. At norm 8 its s is
zero by hypothesis; at norm 10, |s|<=10. For the latter, applying
U and U^(-1) gives |s|<=3q-16, and evenness gives s=q mod 4.
Only s=+-14 need exclusion; the corresponding unit image has
q=8 and |s|=4, contradicting cancellation. The two projected
Pearson inequalities give further necessary conditions that the
normalized dots lie in [-1,1]. These are necessary conditions
only: the survivors need not be realizable in any lattice.

The rational coefficients in delsarte.json define
  F(a,b)=1+sum c_ij [P_i(a)P_j(b)+P_j(a)P_i(b)]/2,
where P_i are normalized Gegenbauer polynomials for S^46. All
c_ij are positive; hence the kernel is positive semidefinite on the
product, by the spherical addition formula and Schur products, and
the sum over all ordered pairs is at least the square of the number
of sites (the constant term is 1). The replay verifies over Q(sqrt(5))
that all 25 surviving labels have negative F. The diagonal is
therefore an upper bound on the cardinality and is exactly
689541416709377/62500000 < 11,040,000.

Run python -B verify.py using Python's standard library. The rational
certificate, not the numerical search, is the proof. The golden
shell hypotheses follow the earlier exact trace-twist reduction in
the 96D golden-trace packet; this result does not assume the
existence of a parent satisfying them.
