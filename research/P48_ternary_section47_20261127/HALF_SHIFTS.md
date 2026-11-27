# A second architecture: half shifts and short dyadic cycles

The 47D benchmark and section notation \(D,v,L\) are those of the [ternary proof](README.md). A nonintegral half shift in a norm-eight section would improve the comparator by \((48/47)^{47/2}-1\approx64.0\%\) if its squared depth were 6. Even a four-cycle could exploit a smaller required depth. These genuinely different gains also fail for the short-normal stock.

**Theorem.** For every \(w\in L\), set \(r=\|w\|^2/2\bmod2\). The raw squared minimum of \(w+2L\) is at most:
\[
\begin{array}{c|cc}
D&r=0&r=1\\ \hline
6,8&20&18\\
10,12,14&20&22\\
16&\text{not asserted}&22.
\end{array}
\]
None of these certificates assumes that \(L+\mathbb Zw/2\) is integral. In particular, the ternary construction's failure cannot be avoided just by going to a short dyadic cycle. For \(D=6,8,10,12,14\), no finite even-order subgroup of at most eight points in \(\mathbb R^{47}/L\) can beat the comparator, and no elementary 2-subgroup of at most 64 points can do so. The assertions concern subgroups and the stated periods only.

## 1. Same Jacobi mechanism, level two

Take the real twist
\[
F_2=\sum_{x\in\Lambda}(-1)^{x\cdot w}q^{\|x\|^2/2}y^{x\cdot v}.
\]
It is even in \(z\) and satisfies the same elliptic law of index \(D/2\). By Poisson, \(F_2|S=\Theta_{\Lambda+w/2,v}\), with integer normal labels. Its exponents are \(q^{r/4}Q^n\), \(Q=q^{1/2}\), \(n\ge0\), because \(\|x+w/2\|^2/2\equiv \|w\|^2/8\bmod\tfrac12\). Thus \(F_2|T=F_2\) and \(F_2|ST^2S^{-1}=(-1)^rF_2\). The eta quotient \(U=\eta(\tau)^8/\eta(2\tau)^4\) has weight 2 and the identical two multipliers \(1,-1\), and leading cusp-zero order \(q^{1/4}\). The two generators and \(-I\) generate \(\Gamma_0(2)\), by the same integer-column argument as at level three (reduce \(|b|\le |u|/2\) then \(|u|\le |b|\)). Therefore \(F_2/U^r\) is weak Jacobi of weight \(24-2r\) and the trivial character. Applying the same evaluation and division proof writes it as \(\sum_j f_jX^jY^{D/2-j}\), with \(f_j\) a scalar modular form of weight \(K=24-2r+2j\).

For these scalar spaces use
\[
A=2E_2(2\tau)-E_2(\tau),\quad
B=\eta(\tau)^{16}/\eta(2\tau)^8,\quad
C=\eta(2\tau)^{16}/\eta(\tau)^8 .
\]
Here \(A\) is holomorphic of weight 2 (the anomalous terms of \(E_2\) cancel), and \(B,C\) are weight 4. At infinity their orders are \(0,0,1\). For \(K=4d+2s\), \(s=0,1\), the monomials \(A^sB^iC^{d-i}\) have distinct initial cusp orders and form a basis by valence, since the index is 3. At cusp zero,
\[
A|S=-\tfrac12 A(Q),\quad B|S=2^4 C(Q),\quad C|S=2^{-8}B(Q),\quad
U|S=-4\,q^{1/4}\prod_{n\ge1}(1-Q^{2n})^{8}(1-Q^n)^{-4}.
\]
The first identity follows directly from the \(E_2\) transformation; the rest follow by the standard eta inversion. These transformations also independently verify holomorphy at zero and the claimed trivial multipliers for \(A,B,C\). The column scalar after extracting \(q^{r/4}\) is \((-1)^j2^{2r-s+4i-8(d-i)}\).

Use the same exact infinity coefficients as before. At cusp zero all coefficients are nonnegative, and \(q^{r/4}Q^n y^b\) is zero if \(2b^2>D(r+2n)\). The coefficients with \(b=0\) and \(2(r+2n)<H\) all vanish under the contrary raw clearance gate \(\min(w+2L)\ge H\). Use \(H=24\), except \(H=22\) for \(D=6,8,r=1\). The standalone [verify_half.js](verify_half.js) computes the exact columns and 11 rational Farkas contradictions, with no external lattice coordinates. Congruence of raw squared norm to \(2r\bmod4\) gives the table.

## 2. Which short enlargements are excluded

A group of \(M\) torus phases with minima \(\mu\) and period \(L\) has squared density ratio \(2M^2(8\mu/47)^{47}/D\). Any even-order subgroup contains a nonzero involution, i.e. \(w/2\) for some \(w\in L\). Its clearance is at most \(22/4\) here. For \(M\le8,D\ge6\), the ratio is at most \(128(44/47)^{47}/6<1\). The bound includes the genuine order-four case regardless of the clearance of the generator.

To treat larger elementary 2-subgroups, \(q(w)=\|w\|^2/2\bmod2\) is a quadratic form on \(L/2L\). On every binary subspace of dimension at least 3 it has a nonzero isotropic vector. Indeed if all nonzeros had value 1, the polar form would take value 1 on every pair of distinct nonzero vectors, contradicting bilinearity on three independent vectors. That phase supplies \(\mu\le20/4=5\). For \(M\le64,D\ge6\) this gives at most \(8192(40/47)^{47}/6<1\). These are exact comparisons replayed by the checker. No general bound is claimed for arbitrary nonadditive real motifs or for a longer normal section. Such changed motifs remain distinct constructive routes.
