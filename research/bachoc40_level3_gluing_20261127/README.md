# Even 3-primary gluing in 40 dimensions: exact density ceiling

**Result.** Let $L$ be the Bachoc–Nebe even 3-modular 40-dimensional lattice, of squared minimum 8 and Gram determinant $3^{20}$. Every even integral lattice $M\supseteq L$ has center density at most $3^{16}/2^{20}$, the **external** Dutour Sikirić et al. 40D benchmark. If equality holds in this family, $M$ has exactly 1,092,000 minimal vectors: the same count as in the external construction. We neither prove equality is attained inside this family nor assert an isometry with that construction. No record on either metric is claimed.

Sources: [Bachoc–Nebe, *Extremal lattices of minimum 8 related to the Mathieu group \(M_{22}\)*](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf), [Nebe's \(L_{40}\) lattice data](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/L40.html), [the authors' exact external 40D benchmark and certificates](https://github.com/wcgbg/sphere-packing/blob/main/README.md). The 3-modular $L$ and the external 40D section are different inputs. The modular-form tools used here are the lattice theta transformation, the eta-quotient criterion, and the valence bound, not a search over lattice classes.

## Exact modular reduction

All norms below are squared. Any even integral $M$ containing $L$ lies in $L^*$. The group $L^*/L$ is 20-dimensional over $\mathbb F_3$, because $L$ is 3-modular. Write $r=[M:L]=3^k$; then $\det(M)=3^{20}/r^2$, so $k\le10$. Put $J=\sqrt3M^*$, a sublattice of $\sqrt3L^*\cong L$ of minimum at least 8. Also $\sqrt3L^*$ is even, so $M$ has level dividing 3. As its determinant is an even power of 3, the weight-20 theta series $\Theta_M(z)=\sum_{v\in M}q^{\|v\|^2/2}$, $q=e^{2\pi iz}$, has trivial character on $\Gamma_0(3)$.

Use the weight-$w$ normalized Fricke involution

$$
(W_w f)(z)=3^{-w/2}(z/i)^{-w}f(-1/(3z)).
$$

Poisson summation gives $W_{20}\Theta_M=r\Theta_J$. For the hexagonal root lattice, $T=\Theta_{A_2}$ and $W_1 T=T$. The weight-3 eta quotients

$$C=\frac{\eta(z)^9}{\eta(3z)^3},\qquad
D=27\frac{\eta(3z)^9}{\eta(z)^3}$$

have the same character $\chi_{-3}$ as $T^3$, and the eta transformation gives $W_3 C=D$, $W_3 D=C$. They are holomorphic at both cusps. Moreover, $T^3=C+D$: their difference has order at least two at infinity (constant and linear coefficients coincide), while the weight-3 valence bound on $\Gamma_0(3)$ permits order at most $3\cdot4/12=1$ for a nonzero form. The seven weight-20 forms $T^{20-3j}D^j$, $0\le j\le6$, have distinct orders at infinity. They span by the same valence bound (order at most $\lfloor20\cdot4/12\rfloor=6$). Thus, with $x=D/T^3$ and $P$ a polynomial of degree at most 6,

$$\Theta_M=T^{20}P(x),\qquad W_{20}\Theta_M=T^{20}P(1-x),\qquad x=27q+O(q^2).$$

## Seven local conditions give the ceiling

If $\min(M)\ge6$, the coefficients at degrees 0 through 2 of the first series are $(1,0,0)$; degrees 0 through 3 of the second are $(r,0,0,0)$. These seven linear conditions determine $P$: a homogeneous difference would be divisible by $x^3(1-x)^4$ but have degree at most 6. Expanding the eta products and $A_2$ theta series to degree four yields

$$[q^3]\Theta_M=1500(r-1),\qquad [q^4]W_{20}\Theta_M=360855-495r=495(729-r).$$

The companion verifier reconstructs these identities over the rationals by solving at $r=0,1$ and taking the affine difference. Thus $r\le729$ from nonnegativity, and $r>1$ gives minimum exactly 6. If $k=0$ the baseline center density is $2^{20}/3^{10}$. If $1\le k\le6$ and $\min(M)=6$, the ratio to $3^{16}/2^{20}$ is $3^{k-6}$. If $\min(M)<6$ it is at most 4 by evenness, and the ratio is at most $2^{20}3^{k-26}\le2^{20}/3^{16}<1$. These cover all $k$. Equality forces $k=6$ and gives $1500(729-1)=1092000$ minimal vectors. The external source reports that same number independently; this conditional count does not identify the lattices.

The equivalent anti-invariant check: the $W_{20}$-anti-invariant basis consists of $T^{20-3j}(D-C)^j$ for $j=1,3,5$. Its degrees 1–3 give row matrix $((-66,-1170,35034),(42,1422,-32034),(150,-7650,110850))$, determinant $-1889568000$. Hence simultaneous minima 8 for $M,J$ force $W_{20}\Theta_M=\Theta_M$, and comparing constants immediately gives $r=1$. The interpolation above also controls minimum 6 and is stronger.

## What could still improve the record

View a nonzero class $u+L$ in $L^*/L$ as a proposed additional phase. The line $0,u,-u$ gives a lattice of index 3; if its minimum stays at 8, its center density is $2^{20}/3^9$, a factor $2^{40}/3^{25}\approx1.2977$ over the incumbent. The exact minimum gate for index $3^k$ is

$$\mu>8\left(\frac{3^{26-k}}{2^{40}}\right)^{1/20}.$$

For $k=1,2$ the right side is between $22/3$ and 8; all norms in $L^*$ are multiples of $2/3$, since $\sqrt3L^*$ is even. Also $L$ itself contains norm-8 vectors, so an index-3 or -9 improvement must have minimum exactly 8. Every nonzero isotropic class has an even integral index-3 extension and therefore has a representative of norm at most 6 by the theorem. Thus an index-3 or -9 improvement must be anisotropic on its nonzero phase differences. More generally, additive subgroups of dimensions 3 through 6 contain a nonzero isotropic class (a degenerate quadratic form has a radical; otherwise diagonalize and use two opposite signs or three equal signs over $\mathbb F_3$), so their density cannot strictly improve either. The argument does not exclude anisotropic lines or planes, larger nonintegral subgroups, or other lattice and nonlinear periodic constructions. No viable anisotropic witness is established.
