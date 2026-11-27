# Exact covering radius for the one-heavy-coordinate shortened Golay stock

Let \(C\) be the even subcode of the perfect punctured binary Golay code of length 23, embedded by its \(0/1\) representatives. Let

\[
L=2D_{23}+C,\qquad Q_\ell(x)=\ell x_0^2+\sum_{i=1}^{22}x_i^2,\qquad 1\le\ell\le2.
\]

The notation \(2D_{23}+C\) means the union of translates by all \(C\) words. It is a lattice: addition differs from XOR by twice the coordinatewise intersection of the two words, and self-orthogonality of \(C\) makes the intersection even. Its covolume in this metric is \(2^{13}\sqrt\ell\), and its squared minimum is 8 (the code has minimum weight 8 and the period has doubled two-coordinate roots supported on light coordinates). **The exact squared covering radius is \((30+\ell)/4\).** Thus the largest center density of two arbitrary translates of this stock, optimized over \(1\le\ell\le2\), is \(1/2\), attained at \(\ell=2\). It does not reach the external record in the main note. Nothing is assumed about rationality or symmetry of the second phase.

## The weighted correction

Fix any target point. Reduce each coordinate modulo 2, then reflect it as needed into \([0,1]\). This introduces only coordinate sign changes and parity shifts: for each word \(c\), its candidate translate after reflection still has binary representative \(c\), plus a parity condition on the multiples of 2. Write the reduced target as \(x_i=1/2+t_i\), put \(s_i=(-1)^{c_i}\), and let \(a_i\) be the weight of coordinate \(i\). The squared residual from the nearest even multiples in this binary phase is

\[
\tfrac14\sum_i a_i+\sum_i a_i t_i^2+\sum_i a_i s_i t_i.
\]

If period parity must be changed, moving the nearest even multiple by 2 at a light coordinate \(i\) costs \(2-4s_i t_i\). Since this correction is nonnegative, it is safe for an upper bound to pay for it even when parity already matches. It suffices to produce a word and a light index \(i\) for which

\[
\sum_j a_j s_jt_j-4s_i t_i \ \le\ -\sum_j a_j t_j^2. \tag{*}
\]

Define \(H=\ell|t_0|\le1,\ S=\sum_{i>0}|t_i|\); always \(\sum a_it_i^2\le (H+S)/2\). A *match* is a binary bit choice contributing \(-a_i|t_i|\); any sign can be specified when \(t_i=0\). The target binary word of all matches has a unique punctured-Golay neighbor \(c_0\) at distance \(e\le3\).

We repeatedly use a small-support estimate. Let \(U\) be the sum of the four largest light \(|t_i|\). Among shortened Golay words, we may prescribe matches on the heavy coordinate and these four lights, and every other bit of the remaining words is still balanced. In fact the dual of \(C\) is the punctured Golay code of minimum 7, so any six prescribed bits are uniform. Fix a light correction coordinate outside the four prescribed ones. The conditional expectation of the left side in (*) is \(-H-U\). Also

\[
\sum a_i t_i^2\ \le\ H/2+U/2+(U/4)(S-U).
\]

It follows that (*) holds if \(S\le2\), or if \(S\le3H\): in the second case even the looser bound \(H/2+U/2+3HU/4\) is at most \(H+U\), using \(U\le2\) and \(H\le1\).

## Finding a correcting word

If \(c_0\) is odd, add a weight-seven Golay word whose support contains all its \(e\) mismatch positions. The heptads form the Steiner design \(S(4,7,23)\), so among heptads through the mismatch set each remaining position participates with probability \((7-e)/(23-e)\). The resulting word is even, and its mismatches are precisely the \(7-e\) new positions. Choose the largest mismatched light coordinate as the correcting index; one always exists. If their absolute values sum to \(V\), the light mismatch cost after the correction is at most \((2-4/(7-e))V\). The heavy mismatch cost, when present, is \(2H\). Averaging over the heptads gives a word with left side in (*) at most

\[
-H-S+\frac{14}{23}H+\frac{10}{23}S=-\frac9{23}H-\frac{13}{23}S.
\]

This is enough when \(S\ge5H/3\). Otherwise \(S<2\) and the small-support estimate applies.

If \(c_0\) is even, the same correction at the largest mismatched light coordinate is enough immediately when there are one or two light mismatches and the heavy bit matches: the corrected light mismatch cost is nonpositive. If the heavy bit mismatches as well, the left side is at most \(-S+H\), which suffices for \(S\ge3H\); use the small-support estimate for \(S<3H\). If there are three light mismatches, the heavy bit matches. The corrected mismatch cost is at most \((2/3)V\le1\), enough when \(S\ge2\); otherwise use the small-support estimate again.

If \(c_0\) has no mismatch, use the smallest light coordinate to correct; the extra cost is at most \(4S/22\), easily absorbed by \((H+S)/2\). The last possibility is a single heavy mismatch. Add a shortened Golay octad containing that position. It leaves precisely seven light mismatches, with each light appearing in a fraction \(7/22\) of such octads. Correct at the largest of them. Averaging gives a left side at most \(-H-S+(2-4/7)(7/22)S=-H-6S/11\), completing (*).

The design facts follow from the extended Golay Steiner system \(S(5,8,24)\): heptads are octads through the deleted coordinate; octads through the heavy coordinate but avoiding the deleted one total 176, and each chosen light lies in 56 of them. For completeness the companion verifier enumerates these finite ingredients from a cyclic generator. The perfectness and uniform-bit properties follow from the punctured minimum 7, linearity and the extended code's self-duality.

## Attainment and density

Take \(w=(1/2,\dots,1/2)+2e_0\). For each \(c\), the uncorrected distances have absolute value \(1/2\) in every coordinate, the correction is required, and its minimum cost is 2 on a light coordinate. The uncorrected norm is \((22+\ell)/4\). This realizes equality in the covering estimate.

The two-translate density ceiling is therefore \(([(30+\ell)/4]/8)^{23/2}/\sqrt{2\ell}\). Its logarithmic derivative has sign \(22\ell-30\), so an endpoint maximizes it. At \(\ell=1\) it is strictly below \(1/2\); at \(\ell=2\) it is exactly \(1/2\). This excludes the fixed-stock two-translate construction for the stated one-coordinate metrics, not arbitrary diagonal metrics, code changes, or additional translates.
