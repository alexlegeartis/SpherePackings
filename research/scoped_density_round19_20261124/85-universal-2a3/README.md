# A universal index-two obstruction for the full family of 2A3 charts

Let \(L\) be any rank-88 even unimodular lattice of squared minimum 8
and let \(K\subset L\) be an embedded scaled root lattice \(2A_3\)
(with Gram 8 on the diagonal and 4 off it). No particular embedding
is assumed. Let \(T=L\cap K^\perp\).
Six of the seven index-two subgroups of \(K^*/K\) and the full phase
group, when used as additive motifs, cannot beat BN85. In fact all
seven covered cases have squared minimum at most \(117/16\).
The seventh index-two subgroup and smaller subgroups are not excluded
by this argument. Nor does it exclude non-additive motifs.

Write \(b_i=x\cdot v_i\) for the three dot labels relative to
the standard star basis \(v_i=2(e_i-e_0)\), and put \(b_0=0\).
The six index-two subgroups have defining parity \(b_i=b_j\bmod2\)
for one unordered pair on \(\{0,1,2,3\}\).
These are actual subgroup conditions because the normal Gram is even.
The action of the coordinate permutations of \(A_3\) is transitive
on the six kernels, so proving one of them for any star basis proves
all six without needing an ambient automorphism.
The remaining index-two group has defining parity
\(b_1+b_2+b_3=0\); these are all seven since the three parities of
the integer dot labels give the full rank-three quotient modulo 2.

For any parent norm-8 vector outside \(K\), all four \(b_i\) differ
by at most 4: test against both signs of each \(K\)-root
\(v_i-v_j\). Each such dot-label is nearest in its phase class, since
subtracting a \(K\)-vector with smaller normal would give total
norm less than 8. The replay checks the finite nearest-label
census and the numerator
\[
 Q(b)=16 b^\mathsf TG_K^{-1}b
      =4\sum b_i^2-(\sum b_i)^2 .
\]
The essential shared constraint is that nearest labels within the
same phase class have *identical shell occupancies*: adding a
\(K\)-vector transports one set bijectively to the other while
keeping horizontal projection and normal norm. Negation identifies
the opposite class too. One may therefore add a dual inequality
over the signed label representatives in a class, not require it
at every representative independently.

The parent shortest-shell size is \(168,498,000\) by
[Sakamoto's theta count](https://arxiv.org/html/2609.15633v2#S6).
The harmonic weighted theta series of degrees 2 and 6
have weights 46 and 50 and each vanishes to order at least 4.
Dividing by \(\Delta^4\) leaves weights \(-2\) and 2, both zero
in level one. Thus the sums of these harmonics on the minimum
shell vanish. The 12 norm-8 vectors in \(K\) are explicit and
removed before applying inequalities to the other labels.

The replay supplies a rational linear combination \(p(b)\) of a
constant and the exact harmonic projections of label monomials
of degrees 2 and 6. It checks that the average of \(p\) over
each nearest signed-label class is at most the indicator of
\(\{b_1+b_2\equiv0\bmod2,\ Q(b)>8\}\). Outside \(K\), summing
using the shared occupancies gives at least \(329747\)
literal norm-8 vectors meeting that indicator, by the harmonic
vanishing and subtraction of the twelve normal vectors.
The nonzero phase with such a vector has \(Q(b)\ge11\),
so its squared separation from the zero phase is at most
\(8-11/16\). Every one of the six subgroup motifs and the full
phase group contains such a forbidden difference.

As before the horizontal period has covolume \(16\), because \(K\)
is primitive by the rank-three Hermite bound and \(L\) is unimodular.
For 128 and 256 phases, the exact comparison
\(M^2\cdot117^{85}<64\cdot128^{85}\) rules out a density
improvement. Everything here concerns specified additive motifs
and not the density of all 85-dimensional packings.

Run python -B verify.py (Python standard library, exact).
