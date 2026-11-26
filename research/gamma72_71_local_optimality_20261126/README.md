# An actual local optimality certificate for the 71-dimensional two-coset packing of Γ72

No density record is claimed. Start from the explicit norm-10 normal in [the original 71D certificate](../novel_routes_20261113/higher_parent/gamma72_verified.json): the center density is (2(79/40)^{71/2}/\sqrt{10}=19,663,555,319.937\ldots). Its published background is [Nebe's minimum-8 theorem](https://arxiv.org/abs/1008.2862) and the [Conway–Sloane antipode construction](https://arxiv.org/abs/math/0207182). This result concerns the actual contact geometry, **not a normal-diameter estimate**. In the usual two-translate parameterization, after dividing out a common similarity, changing the entire period Gram and moving either or both sites cannot improve it within the explicit neighborhood below. Larger-period perturbations, distant metrics, and other configurations are not covered.

## Contacts and the perturbation to be attempted

Let an even unimodular lattice of dimension 72 have minimum 8, and fix a lattice vector v of squared norm 10. It is primitive. Let P=the perpendicular section Γ ∩ v^⊥. Take the two projections of the integer sections at heights x·v=0 and 1. The first is P, the other is c+P; both heights occur by primitivity and unimodularity. P has covolume √10. For a shell point x (squared norm 8) write t=x·v and y=x−tv/10. Then |t|≤5 by the minimum applied to x±v. Every shell point at t=1 gives a true cross-contact of squared length μ=79/10. No shorter difference is possible by lifting in the two sections. At least one such contact exists by the positive second moment below. Within either section distances squared are ≥8.

Work with an arbitrary traceless symmetric H on the 71D perpendicular space and a linear form l(y)=a·y. Write Q(y)=yᵀHy, h=||H||F, u=||a||, and N²=h²+u². We prove for some actual contact y with t=1 that
`Q(y)+l(y) < −N/100000` if N>0.
Thus even correlated shape and motif changes lose a contact; single-motif testing would not suffice.

## Moment certificate with no enumeration of the 72D shell

Use expectation on the entire normalized shortest shell, not on a conditional probability space. It is an 11-design: in degrees 2,4,6,8,10 the harmonic theta form has weight 36+j and vanishing order ≥4, so division by Δ⁴ yields negative weight j−12. Odd moments vanish by antipodality. This is the classical design property; see [Bannai et al.](https://doi.org/10.17323/1609-4514-2006-6-2-225-264). For any homogeneous tangent polynomial W of degree d and polynomial p(X)=Σ p_k Xᵏ with total degree 2deg(p)+d≤10, spherical integration gives
`E[p(t²)W(y)] = E[W(y)] R_d(p)`,
`R_d(p)=Σ p_k 80ᵏ (2k−1)!! / Π_{i=0}^{k−1}(72+d+2i)`.
This follows by independent Gaussian coordinates along v and its perpendicular, then dividing by their joint radial moment. For traceless H: E[Q]=0, E[Q²]=128h²/(72·74), E[l²]=8u²/72.

Write [j] for the indicator of |t|=j. Let p_{j,S}(X)=Π_{s∈S}(X−s²)/(j²−s²). The following rational inequalities hold pointwise at X=0,1,4,9,16,25; exact arithmetic is reproduced in verify.py.

* p5=p_{5,{0,1,2,3,4}}=[5], with E[5]=142/77727195. P5⁺=p_{5,{0,1,2,3}} majorizes [5], R₂(P5⁺)=83/2878785.
* P₂=p_{1,{0,2,3,4}}=[1]−210[5], and R₂(P₂)=39952/82251. P₄=p_{1,{0,2,5}}≤[1], with R₄(P₄)=850/2223.

Consequently E[[1]l²]≥C u², E[[1]Q²]≥B h², where C=39952/740259 and B=6800/740259. Also E[[1]Q]=210 E[[5]Q]. To bound its bias, M=E[[5] yyᵀ] is positive semidefinite with trace T=(11/2) E[5] (since |y|²=8−25/10), and largest eigenvalue at most U=(8/72)83/2878785 by the majorant. Since 3U<T<4U, convexity gives tr(M²)≤3U²+(T−3U)². Tracelessness of H thus gives
`|E[[1]Q]| ≤ 210 sqrt(3U²+(T−3U)²−T²/71) h < A h,`
where A=571/500000.

Orient all heights −1 back to +1 by replacing their y by −y. The oriented linear and cubic tangent averages vanish. For a linear tangent polynomial use an odd polynomial in t equal to 1 at t=1 and 0 at t=2,…,5, e.g. t times the even Lagrange selector on their squares, of degree 9. Its product with the linear form has degree 10 and zero spherical expectation. For a cubic, use the degree-7 selector zero at t=2,3,4: the leftover fifth shell contributes exactly zero, since x↦v−x pairs its height-five tangents y and −y; the height−5 analogue is −v−x. Spherical expectation still vanishes. In particular the cross-term Ql vanishes in the oriented second moment. For F=Q+l on the oriented contacts and E₁ the unnormalized mass with |t|=1, we obtain |E₁F|≤Ah and E₁F²≥Bh²+Cu².

Put S=1/10000. Exact arithmetic gives B−Aμ>S, C>S and
`4(B−Aμ−S)(C−S) > A² μ.`
Hence `Bh²+Cu² > Ah(μh+√μ u)+SN²` for N>0. On contacts F≤M₀=μh+√μ u<9N. If F≥−ηN everywhere for η=1/100000, then F²≤M₀F+ηN(M₀+ηN) pointwise, implying E₁F²≤M₀Ah+(9+η)ηN², in contradiction since (9+η)η<S. This proves the claimed quantitative contact witness.

## Passage from the linear witness to the actual nonlinear packing

Model an arbitrary nearby full period and two sites by a positive definite Gram G relative to the original metric, plus an independent relative translation b of the second site in the pulled-back coordinates; discard common translation, rotation and overall scale. Normalize det G=1 and write G=exp H; then tr H=0. Use a=2b and N²=||H||F²+4||b||². At a cross-contact y the new squared length is (y+b)ᵀexp(H)(y+b). After its first-order term Q+l, the error is at most
`exp(N) (μ/2+√μ+1/4) N² < 8N²` for 0<N≤10⁻⁶.
The offending contact is strictly shorter than μ. The period covolume was fixed, hence the actual density strictly decreases. This is strict local optimality modulo similarities for the fixed two-site periodic description, with an explicit (very conservative) certified radius. It does not prove full supercell jamming.

We also carried out a genuinely different constructive test: split both sites in an index-two supercell with character χ of P. For a fixed Gram an alternating velocity can survive *only if* one of the two character halves of the oriented contact differences fails to span the 71D perpendicular space: for contact y=c+p the two necessary bar equations are y·(b₁−b₀)=0 in χ(p)=+1 and y·(b₁+b₀)=0 in χ(p)=−1, from the two base vertices and conjugate signs. Conversely a missing span supplies a first-order flex of that mode, not automatically a packing improvement. This is a precise search criterion, not a claimed supercell flex.

For the actual pinned normal in Nebe's catalogue basis (v=e44+e51, zero-based), we tried the 72 characters “parity of the ith parent coordinate” restricted to P. Three rounds of applying the six published, exactly checked isometries to norm-eight basis vectors give 7,796 distinct shell points and 2,721 distinct oriented height-one contacts. For each coordinate, both halves of the oriented contact set have rank 71 already modulo 101, on their exact integer rows 10x−v. Therefore no one of these 72 parity choices admits a nonzero alternating velocity. This is a concrete attempt to escape the local theorem by changing periods, not an exclusion of all index-two supercells. Run `OPENBLAS_NUM_THREADS=1 python index2.py` with NumPy; an explicit path to the pinned existing catalogue file is optional. Exact source SHA-256: 92daeaf0d93672e82ca5e5bc4479e5775daf59029aa4e983fb464e36ef002044.

**Strict local consequence for each tested doubled cell.** Allow the period metric and all four site positions to vary in each such fixed index-two description; normalize the determinant. If every active contact has nonnegative first derivative, average the two copies of every edge over the two quotient vertices. This gives exactly the primitive two-site inequalities for H and the mean relative site velocity. The strict witness above forces H=0 and those two mean velocities equal. Subtract their common translation. Each averaged edge now has derivative zero, so both of its nonnegative constituent derivatives are zero. Only the alternating velocities remain, and the two full ranks just verified force both to vanish. There is therefore no nonzero first-order permissible tangent for each of these tested four-site descriptions. Compactness of its unit tangent sphere and the finite contact set yield strict local optimality modulo similarities in some neighborhood of that representation. No uniform radius for these four-site neighborhoods, and no statement for arbitrary supercells, is asserted.
