# Independent 15D route: the low-height Barnes–Wall cut gate

Let \(B=(C+2D_{16})/\sqrt2\), where \(C=\mathrm{RM}(1,4)\) in
its 16-bit integral representatives, and \(D_{16}\) is the
integer even-coordinate-sum lattice. This is the usual
[Construction B/low-dimensional Barnes–Wall model](https://neilsloane.com/doc/Me14.pdf)
with squared minimum 4 and covolume 16. Its norm-four shell
consists of the numerator vectors \(2(\pm e_i\pm e_j)\) and
the \(\pm1\) vectors on affine binary hyperplanes, with an even
number of negative entries. Both types are immediate from the
Construction B congruence and attain norm eight before
dividing by \(\sqrt2\).

Consider a primitive dual normal \(u\in B^*\), let \(n=u^2\),
and set the period \(L=B\cap u^\perp\). Choose two phases of
\(L\) from two full parallel cut heights whose pairing difference
is an integer \(r\). Then the covolume is \(16\sqrt n\), since
the pairing map is onto \(\mathbb Z\), so height spacing is
\(1/\sqrt n\). At minimum squared separation at most 4, a
density gain over \(1/(16\sqrt2)\) needs \(n<8\). The following
gate covers all \(1\le |r|\le\sqrt{2n}\) in that range.

**Claim.** For every dual normal of squared length \(n<8\)
and every integer \(1\le |r|\le\sqrt{2n}\), some norm-four
vector \(x\in B\) has \(\langle u,x\rangle=r\).

To see this, duality against \(2D_{16}\) and \(C\) in raw
coordinates gives \(u=y/\sqrt8\), with \(y\) all even or all
odd, and \(y\cdot c=0\bmod4\) for all \(c\in C\). Hence
\(n=y^2/8\) is integral: for even \(y=2v\), the parity
\(v\bmod2\) orthogonal to \(C\) has even weight; for odd \(y\),
the 16 squares are each 1 modulo 8. Each nonzero normal has
\(n\ge2\) from these same constraints: even \(v\) of
squared length <4 and nonempty odd support would contradict
orthogonality to affine functions (the dual binary minimum
weight is four), and nonzero all-even \(v\) has length squared
at least four. Odd \(y\) always has length squared at least 16.

The first norm-four shell type realizes all \(|v_i\pm v_j|\)
if \(y=2v\); for odd \(y\) it realizes
\(|y_i\pm y_j|/2\). The finite count of magnitude patterns in
the accompanying script shows the only missing low \(r\)'s:
single, two, or three nonzero \(|v_i|=2\);
support sizes 10,12,14 with \(|v_i|=1\) on the support;
and the all-absolute-one odd \(y\).

For one or two magnitudes 2, an affine hyperplane can isolate
one nonzero magnitude. For three it can contain all three;
their signs may be chosen to produce both 1 and 3, and one
unused position on the hyperplane can repair sign parity.
For support sizes 10,12,14 with magnitudes 1, use an affine
hyperplane meeting the support in six points. All intersection
counts are even because the support is in \(C^\perp\).
For support 10, if count 6 never occurs then count 4 does not
occur by complements; the only possible counts are 2 and 8.
All 15 nontrivial binary Fourier coefficients of the support
would have magnitude 6, contradicting Parseval
\(10^2+15\cdot36>16\cdot10\).
For support 12, absence of count 6 leaves only 4 or 8,
again contradicting \(12^2+15\cdot16>16\cdot12\).
For support 14 the two possible counts are 6 and 8.
Align the signs of the six active coordinates on the
selected hyperplane and use an unused item to fix even
sign parity, giving pairing 3. For all-absolute-one odd
\(y\), its negative sign mask is in \(C^\perp\) by
the congruence \(y\cdot c=0\bmod4\); the signs can
therefore be fully aligned on any hyperplane, giving 2.

Any such vector produces an actual horizontal cross difference
of squared norm \(4-r^2/n\); there is no nearest-lift assumption.
The ratio of the two-phase density to the benchmark is
at most
\[
 2\sqrt{2/n}(1-r^2/(4n))^{15/2}
 \le \sqrt2(15/16)^{15/2}<1.
\]
For the second inequality use \(r\ge1\); the one-variable
expression has its maximum at \(n=4\).
The strict last inequality follows from
\((16/15)^{15}>1+15/15=2\).
This only rules out the low-height two-phase mechanism:
higher pairing differences, more phases and new normal
metrics are not covered. No new packing or density is asserted.
