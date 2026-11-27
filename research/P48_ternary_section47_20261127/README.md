# Ternary lifts of short primitive sections of extremal 48D lattices

Sphere147, 27 November 2026. **No new packing-density record.** The comparator is the [Chen–Hu–Li–Wang–Wu 47D antipode](https://arxiv.org/html/2505.02394v1#S3.SSx5), center density \(\delta_*=47^{47/2}/2^{118}\), corresponding to occupied density \(2.26696032293109\cdot10^{-8}\) in the [project table](https://github.com/alexlegeartis/SpherePackings/blob/main/LOWER_BOUNDS_1_96.md). The construction tested here is not the paper's: a **nondual** ternary shift in one short-normal section, repeated to give a cyclic three-site motif. In the norm-eight section, the previously identified target \(\min(w+3L)\ge52\), if attainable, would give a 1.34% density gain. We prove that target impossible, and also exclude the four other short normal lengths.

**Theorem.** Let \(\Lambda\) be an even unimodular lattice in dimension 48 of squared minimum 6. Take \(v\in\Lambda\) with squared norm \(D\in\{6,8,10,12,14\}\) and \(L=\Lambda\cap v^\perp\). For any \(w\in L\) put \(r=\|w\|^2/2\bmod3\), in \(\{0,1,2\}\). The squared minimum of \(w+3L\) is bounded by
\[
\begin{array}{c|ccc}
D &r=0&r=1&r=2\\ \hline
6&48&44&46\\
8,10&42&44&46\\
12&48&44&46\\
14&48&50&46
\end{array}
\]
(This allows the zero when the phase is trivial.) Consequently no packing of three distinct translates \(L\cup(L+w/3)\cup(L-w/3)\) in this family improves the stated comparator. Norm sixteen is additionally excluded for residues \(r=1,2\), with bounds 50 and 52 respectively; no uniform statement is made here about its residue \(r=0\) or larger normals. Arbitrary continuous shifts, noncyclic three-site motifs, partial fibers, and changes to the period are outside the theorem.

## 1. Geometry and gain threshold

Such \(v\) are primitive: a proper divisor has squared norm at most \(D/4<6\). Unimodularity then makes \(x\mapsto x\cdot v\) onto \(\mathbb Z\), since otherwise \(v\) would be divisible in \(\Lambda^*=\Lambda\). Slices of height \(1/\sqrt D\) show \(\operatorname{covol}(L)=\sqrt D\). For nontrivial \(w\bmod3L\), the three sites form the lattice \(L+\mathbb Z w/3\) with covolume \(\sqrt D/3\). Its squared minimum \(\mu\) is at most \(\min(w+3L)/9\). The squared center-density ratio is exactly
\[
R^2=\frac{18}{D}\left(\frac{8\mu}{47}\right)^{47}.
\]
At raw norm 52 and \(D=8\), the ratio at \(\mu=52/9\) is \(R^2=9\cdot416^{47}/(4\cdot423^{47})>1\). The theorem's upper bounds give \(R<1\) in every stated case, by integer arithmetic. Trivial phases are irrelevant to the three-distinct-site search. 

## 2. A Jacobi twist that preserves the normal coordinate

Use \(q=e^{2\pi i\tau},\,y=e^{2\pi iz}\), \(m=D/2\), \(a=w/3\), and \(\zeta=e^{2\pi i/3}\). Define the even-in-\(z\) twisted Jacobi theta series
\[
F(\tau,z)=\sum_{x\in\Lambda}\frac{\zeta^{x\cdot w}+\zeta^{-x\cdot w}}2 q^{\|x\|^2/2}y^{x\cdot v}.
\]
It has weight 24 and index \(m\) under the usual elliptic transformations: adding an integral multiple of \(v\) to \(x\) leaves the twist unchanged since \(w\cdot v=0\). Use the standard Jacobi slash operator, including its exponential factor. Evenness of \(\Lambda\) gives \(F|T=F\). Poisson summation gives
\[
 F|S=\tfrac12(\Theta_{\Lambda+a,v}+\Theta_{\Lambda-a,v}).
\]
Both shifted theta series use the same \(y^{x\cdot v}\) for a shifted vector \(x\). For completeness, completing the quadratic Gaussian against the original linear twist turns \(\Lambda\) into its dual translated by \(\pm a\); the Gaussian gives \((-i)^{24}=1\) after the slash weight, and the Jacobi exponential is canceled by the index \(m=\|v\|^2/2\). The real symmetrization makes the sign immaterial. There is no cross-term shift in the \(y\) label because \(a\cdot v=0\).

The cusp-zero exponents of \(F|S\) are \(q^{r/9}Q^n\), \(Q=q^{1/3}\), \(n\ge0\), since \(\|x\pm w/3\|^2/2\equiv \|w\|^2/18\bmod\tfrac13\). Thus \(F|ST^3S^{-1}=\zeta^r F\). Let \(U=\eta(\tau)^3/\eta(3\tau)\). Dedekind's transformation yields \(U|T=U\) and
\[
U|S=(-i)\sqrt3\,\frac{\eta(\tau)^3}{\eta(\tau/3)},\qquad U|ST^3S^{-1}=\zeta U .
\]
The first column \((u,3b)\) of a member of \(\Gamma_0(3)\) can be reduced by these two generator types: a power of \(ST^3S^{-1}\) reduces \(|b|\le|u|/2\), then a power of \(T\) reduces \(|u|\le3|b|/2\) if \(b\ne0\); here \(u\not\equiv0\bmod3\). With \(-I\), this proves the generation assertion. Hence \(F/U^r\) is a weak Jacobi form on \(\Gamma_0(3)\) of weight \(24-r\), index \(m\), character \(\chi_{-3}^r\), even in \(z\). It is holomorphic on the upper half-plane; at both cusps it has nonnegative exponents in their local parameter, as \(U\) has constant term one at infinity and leading order \(q^{1/9}\) under \(S\). There are just the two cusps of \(\Gamma_0(3)\).

## 3. Small explicit basis, with both cusp transformations

Let \(X=\phi_{-2,1}=-\theta_1(\tau,z)^2/\eta(\tau)^6\), \(Y=\phi_{0,1}=4\sum_{k=2}^4\theta_k(\tau,z)^2/\theta_k(\tau,0)^2\) be the standard weak Jacobi forms (see [Eichler–Zagier](https://doi.org/10.1007/978-1-4684-9162-3) and the explicit [theta formulas](https://people.mpim-bonn.mpg.de/zagier/files/arxiv/1208.4074/1208.4074v1.pdf)). Then
\[
 F/U^r=\sum_{j=0}^m f_j X^jY^{m-j},\qquad f_j\in M_{24-r+2j}(\Gamma_0(3),\chi_{-3}^r).
\]
This decomposition applies in this congruence case directly: evaluate at \(z=0\), where \(Y=12\) and \(X=0\), subtract the \(j=0\) term, and divide by \(X\). The numerator is even in \(z\) so has a double zero there, and elliptic transformations transport it to every zero of \(X\). The quotient stays weak at both cusps: their leading Laurent coefficient of \(X\) is \(y+y^{-1}-2\), nonzero. Repeating reduces to index zero, which has no elliptic variable. This avoids assuming a level-one ring theorem automatically holds at level three.

Write \(A=\Theta_{A_2}\), \(B=U^3\), \(C=\eta(3\tau)^9/\eta(\tau)^3\). All have character \(\chi_{-3}\), weights \(1,3,3\); this follows from the eta formulas and the Poisson theta formula for \(A_2^*\simeq A_2/\sqrt3\). At infinity their orders are \(0,0,1\). For \(K=24-r+2j=3d+s\) with \(0\le s<3\), a basis of the indicated scalar space is
\[
 A^s B^i C^{d-i},\qquad 0\le i\le d.
\]
The weights and characters are correct since \(s+d\equiv r\bmod2\); their distinct leading orders at infinity prove independence, and the valence bound at infinity \(K[\mathrm{SL}_2(\mathbb Z):\Gamma_0(3)]/12=K/3\) bounds the dimension by \(d+1\).

Under \(S\), in the variable \(Q\), the three scalar factors become \((-i)3^{-1/2}A(Q)\), \((-i)^3 3^{3/2}C(Q)\), and \((-i)^3 3^{-9/2}B(Q)\). The weak generators \(X,Y\) are invariant with their own slashes, so they use their ordinary series at \(q=Q^3\). After extracting \(q^{r/9}\) from \(U^r|S\), its remaining scalar product in \(Q\) is \(\prod_{n\ge1}(1-Q^{3n})^{3r}(1-Q^n)^{-r}\), multiplied by \((-i)^r3^{r/2}\). Thus a column \((j,i)\) has cusp-zero scalar \((-1)^j3^{(r-s+3i-9(d-i))/2}\); the power is integral. At infinity its factors are the ordinary \(U^r,A^s,B^i,C^{d-i},X^j,Y^{m-j}\). These expansions give a finite exact linear model.

## 4. Rational contradiction to the missing horizontal coefficients

At infinity, \([q^0]F=1\) and \([q]F=[q^2]F=0\). In shell \(q^n\), if \(x\ne\pm v\) then \(\|x\pm v\|^2\ge6\) forces \(2|x\cdot v|\le D+2n-6\). Beyond this range only the specified representative \(x=v\) for \(y^D q^{D/2}\) contributes among nonnegative labels; its twist is one. We use only the first four powers \(q^1\) through \(q^4\) (through \(q^6\) in the indicated norm-16 residue-one certificate).

At cusp zero all coefficients are nonnegative, and \(q^{r/9}Q^n y^b\) is zero if \(9b^2>2D(r+3n)\), simply by Cauchy–Schwarz. Crucially, the **\(b=0\) terms** count \((L+w/3)\) and its negative and nothing with nonzero normal height. If \(w+3L\) had raw squared minimum at least \(H\), then the further coefficients with \(b=0\) and \(2(r+3n)<H\) would all vanish. For each table entry, use the contrary gate equal to that entry plus six, or an equivalent smaller even cutoff immediately above the entry (there is only one possible norm residue modulo six). The literal case list in the replay gives these cutoffs.

[verify.js](verify.js) evaluates these columns in exact integer \(q,y\) arithmetic, imposes the stated equalities, and verifies a finite rational Farkas relation for each case. Specifically, each selected nonnegative cusp term gives a linear inequality in the coefficient variables. Eliminate the equalities over \(\mathbb Q\). A positive linear combination of the selected inequalities cancels every remaining variable and has a strictly negative constant on the right of \(0\le\text{constant}\), a contradiction. The support consists of listed literal ordinary cusp coefficient labels; all multipliers and RHS checks are reconstructed with BigInt fractions, not solver tolerances. The largest cusp-zero local series index used is 12. The resulting strict raw bound reduces to the stated residue-specific table because norms of \(w+3L\) are even and congruent to \(2r\bmod6\).

## 5. Consequences for genuinely nonlinear ternary motifs

This excludes more than cyclic triples. In \((\tfrac13 L)/L\) any two different phases differ by one of the bounded shifts, so the compatibility upper raw squared norms are \(48,46,46,48,50\) for normal squared lengths \(D=6,8,10,12,14\), respectively. For an arbitrary \(M\)-site motif with raw squared clearance at most \(d\), its squared density ratio is at most
\[
\frac{2M^2}{D}\left(\frac{8d}{423}\right)^{47}.
\]
The first motif sizes that could conceivably improve the comparator are respectively **17, 53, 60, 24, 10**; all smaller nonlinear motifs lose. The arithmetic compares both consecutive integers exactly. The four-site almost-regular candidate over \(D=12\) or \(14\) at raw clearance 52, which would have had positive gain, is therefore ruled out as well. This does not assert that any of the listed sizes can be achieved. In the \(D=14\) route with raw clearance 50 the table requires every pair difference to have ternary quadratic value 1, giving a precise finite target for a genuinely larger motif rather than another three-cycle.

An independent [level-two replay](HALF_SHIFTS.md) also rules out the short even-order subgroup route, without assuming shifts are integral.

Replay: `node verify.js`. No external lattice coordinates, shell sampling, floating optimization, or unbounded CVP is used.
