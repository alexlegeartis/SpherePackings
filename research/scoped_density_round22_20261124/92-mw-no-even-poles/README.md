# MW92: no pole-centered even numerator for depth four and norm 26

For a depth-four norm-26 section use the rational pole form from
the accompanying sixth-jet update. Its unique finite pole \(a\)
belongs to \(\mathbb F_{47^2}\). Write \(z=t-a\), \(s=z^2\); after
clearing denominators the numerators \(P(z),Q(z)\) have degrees
26 and 39 and satisfy
\[
Q^2-P^3=z^6(z^{47}-z+a^{47}-a),
\qquad P(0)Q(0)\ne0 .
\]
**The numerator \(P(z)\) cannot be even.** This treats parity
around the pole \(a\), not around a general center, and applies
to both split and nonsplit finite poles. It does not imply that
any additional sixth-jet orbit is excluded.

Suppose \(P(z)=P_0(s)\), of degree 13, and decompose
\(Q=A(s)+zB(s)\). The exact odd part forces
\[
2 A(s)B(s)=s^3(s^{23}-1).
\]
Since \(\deg Q=39\), \(\deg B=19\) and \(\deg A=7\).
The nonzero \(Q(0)\) ensures \(A(0)\ne0\), so the factor \(s^3\)
lies entirely in \(B\). Write \(A=\alpha A_0\) and
\(B=(2\alpha)^{-1}B_0\), with \(A_0,B_0\) monic. The 23 distinct
roots of \(s^{23}-1\) all lie in \(\mathbb F_{47}\); \(A_0\)
takes exactly seven roots, with the other 16 forming
\(s^{-3}B_0\).

Let the superscript \({}^*\) denote reversal at each polynomial's
stated degree, and make \(P_0\) monic by dividing by its leading
coefficient. Calling the monic result \(H\), the even part
implies that its reversal must satisfy
\[
(H^*)^3=(B_0^*)^2+O(u^{25})
\]
at \(u=0\). Indeed \(sB_0^2\) has degree 39, while the \(A_0^2\)
term has degree 14 and the pole-dependent even term has degree
3. Thus the first 25 coefficients of \(H^*\) are those of the
unique formal series \(C\) with \(C(0)=1\) and
\(C^3=(B_0^*)^2\) to that order. If \(B_i\) and \(C_i\) denote
the coefficients of \(B_0^*\) and \(C\), respectively,
differentiating gives \(3 B_0^* C'=2 (B_0^*)' C\), hence
\[
 C_n=(3n)^{-1}\sum_{i=1}^{\min(16,n)}
             (5i-3n) B_i C_{n-i},\qquad 1\le n\le24 .
\]
These equations are over \(\mathbb F_{47}\) because \(B_0^*\)
has coefficients there. Since \(\deg H^*=13\), every coefficient
\(C_{14},\ldots,C_{24}\) must vanish.

The finite replay enumerates all \(\binom{23}{7}=245,157\)
root assignments. Their first nonzero disallowed coefficient
is at degree 14 for 239,476 assignments, degree 15 for 5,612,
and degree 16 for 69. None survives. The test is independent
of the unknown scalar \(\alpha\) and of the finite pole.

Compile and replay:

    c++ -std=c++20 -O3 evenpole.cpp -o evenpole
    ./evenpole

Use assertions enabled (no `-DNDEBUG`). The finite verifier needs
no number-field library or numerical lattice reduction.

### No nontrivial rotation about the pole

Even over the algebraic closure a putative multiplicative symmetry
of P(z) about z=0 has order dividing 26, by its leading coefficient
and nonzero constant. Orders 2 and 26 are excluded above. If the
order were 13, the monic reversal P* would have exponents only
0,13,26. The normalized reversed equation then has the form
(Q*)²=(P*)³+βu^25+O(u^71), with β nonzero. Taking the unique
formal square root with constant 1, the unperturbed term has
exponents in 13Z; the linear term has exponents in 25+13Z.
The quadratic term has first coefficient −β²/8 at u^50, before
any cubic or residual term. This contradicts degree 39 of Q*.
