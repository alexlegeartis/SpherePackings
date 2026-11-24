# Five-torsion lattice extensions of Kappa11 cannot increase density

Let \(K\) be the integral [Kappa11 of Nebe–Sloane's catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KAPPA11.html), with Gram determinant 972 and minimum squared norm 4. Every lattice \(K\subseteq L\subseteq (1/5)K\) has center density at most \(1/\sqrt{972}\), attained by \(K\) itself. This does **not** constrain arbitrary 11D lattices, arbitrary periodic motifs in this torus, other denominators, or shape changes. It does not change a lower or an unrestricted upper bound.

## Finite formulation

Put \(C=5L/5K\subset K/5K\cong\mathbf F_5^{11}\), an additive subgroup of order \(5^r\). For a residue \(x\), let \(m(x)\) be the least squared norm of an integral lift. For a nonzero residue define
\[
m_*(x)=\min(m(x),m(2x)).
\]
The minus sign makes no difference, so this is the least lift norm along the full nonzero projective line. If the shortest squared norm of \(L\) is \(d/25\), every nonzero \(x\in C\) has \(m_*(x)\ge d\); also \(d\le100\) by the period. All possible \(d\) are even by integrality and the even diagonal of the Gram. Covolume is \(\sqrt{972}/5^r\); beating \(K\) would require
\[
 5^{2r}d^{11} > 100^{11}.
\]
Ranks 0 and 11 give \(K\) and \((1/5)K\), which have the same density. The replay computes \(\max_{x\ne0}m_*(x)=48\). That already excludes ranks 1 and 2, whose improvement thresholds are 76 and 56.

## Integer Fourier certificates at the other ranks

In the Fourier bound we simultaneously use both consequences of scalar symmetry: the admissible set uses \(m_*\), and Fourier characters are averaged over full nonzero scalar lines. Four checked integral isometries (transpose action on characters) and multiplication by 2 generate 714 point orbits and 714 character orbits, including their identity orbits. The character sum of a single orbit \(O\) evaluated at a given point is exactly \((5N_0-|O|)/4\), where \(N_0\) is the number of its characters annihilating the point. This follows because each of the four nonzero pairing values occurs equally often, and their fifth roots of unity sum to \(-1\). This collapses the earlier, larger real-cyclotomic system to an integer one without discarding any necessary conditions on \(C\).

For each row below, the supplied nonnegative integer weights on nontrivial character orbits yield a positive definite \(F\) with no constant coefficient, such that \(F(x)\le-P<0\) at every nonzero allowed difference with \(m_*(x)\ge d\). Sum \(F\) over all ordered motif differences. The Fourier expression is nonnegative, whereas the diagonal contributes \(|C|F(0)\) and the remaining terms at most \(-P|C|(|C|-1)\). Therefore \(|C|\le1+F(0)/P\). The replay derives \(P\) and the following caps by exact integer arithmetic; numerical optimization is not part of verification.

| Rank | First improving d | Cap on motif size | Required size |
|---:|---:|---:|---:|
|3|42|60|125|
|4|32|410|625|
|5|24|2927|3125|
|6|18|14639|15625|
|7|14|48706|78125|
|8|10|297868|390625|
|9|8|805394|1953125|
|10|6|2682882|9765625|

Ranks 9 and 10 also admit a separate short-norm certificate: a quotient map to \(\mathbf F_5^{11-r}\) cannot avoid zero on the short nonzero residue directions. The search represents its column labels in order of first appearance of independent labels, setting each newly independent column to the next unit vector by a change of codomain basis. The last nonzero coordinate of a short direction forbids precisely one label, once the preceding labels are fixed. Rank 9 (\(d=8\)) prunes after nine input coordinates; rank 10 (\(d=6\)) prunes after four. These exact branches are supplementary to the Fourier proof.

## Replay

Run in this directory with the C++17 standard library and Python 3 standard library:

    python check_inputs.py
    g++ -O3 -std=c++17 short.cpp -o short
    ./short 4 0 0 < K11.ldl
    ./short 6 1 11 < K11.ldl
    ./short 8 2 11 < K11.ldl
    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify

The replay traverses all character and point classes; it regenerates pairings from the lattice isometries rather than trusting a stored Fourier matrix. It uses exact nearest-first LDL searches on orbit representative lifts and their doubles rather than a complete archived CVP table. There are no numerical solver requirements. The original metric is independently checked by reconstructing the Gram from the triangular identity and by checking both the literal norm-4 vector and the exhaustive exact short-vector search. The second and third invocations of short are independent checks for ranks 10 and 9.
