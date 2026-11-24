# Ternary stacking obstructions and the paired depth of fifth-torsion

These results concern restricted layer constructions and simultaneous depths, not the unrestricted packing or covering optimum. Let L be an extremal even unimodular lattice of rank n=24s and squared minimum m=2s+2, with D=2s. The main cases are n=48,72. Assertions for n=96 are conditional on an extremal even unimodular L of that dimension.

## Third-torsion data used here

For r=v²/2 mod 3, the cap on d=dist(v/3,L)², the nearest count if equality holds, and the published proof of the cap in ranks 48 and 72 are as in the prior [third-torsion certificate](../../order_three_tori_48_72_20261123/README.md). The exact replay here also handles 24 and conditional 96. The additional least-shell coefficient of the twisted theta series, when a phase equals its residue's depth cap, is as follows:

| n | depth caps (r=0,1,2) | coef of q^(s+1) at caps (r=0,1,2) |
|---:|:---|:---|
| 24 | 2, 14/9, 16/9 | −270, −27, −108 |
| 48 | 4, 32/9, 34/9 | −19,512, −5,904, −11,250 |
| 72 | 6, 50/9, 52/9 | −1,278,423, −480,411, −802,872 |
| 96 | 8, 68/9, 70/9 | −82,872,720, −34,230,924, −54,117,315 |

For self-contained replay, the eta basis underlying the prior proof is:
write n/2−r=3J+a, 0≤a<3, put A=theta(A2), U=eta(τ)^3/eta(3τ), V=eta(3τ)^9/eta(τ)^3. The twisted form lies in the span of A^a U^(3j+r)V^(J−j), 0≤j≤J; their local orders at 0 are j+r/3, and at infinity J−j. Indeed T fixes the twist, ST³S⁻¹ acts by exp(2πir/3), and U has the r=1 multiplier; dividing by U^r leaves an ordinary form with finite character and weight n/2−r. The A,U³,V products are holomorphic at both cusps with matching character; valence at index 4 gives the dimension bound. Initial coefficients (1,0,…,0) through q^s determine the surviving j≥J−s basis coefficients uniquely at the cap. The replay expands one extra coefficient to obtain the table.

If the r=0 cap is attained, harmonic theta forms of degree h=1,2 vanish by valence: at infinity order at least s+1, at zero at least 3s, larger in sum than (12s+h)/3. The possible odd sign can be removed by squaring. Thus each nonempty shifted shell is a spherical 2-design (no antipodality needed). Also, if every nonzero element of a linear ternary phase subspace has depth D, Fourier inversion of the coefficients at q^(s+1) shows 3^dim divides N−a0, where N is the parent least-shell size and a0 is the r=0 coefficient in the table. This bounds the dimension by 9,11,13,11 respectively. Specifically at any nontrivial evaluation pattern, exactly (N−a0)/3^dim original shell vectors realize that pattern.

## New ternary stacking bounds

Consider full horizontal translates of L, with all differences of horizontal phases in (1/3)L/L. Let the packing's squared minimum be μ>0 and vertical dimension k; μ≤D+2 because each horizontal layer contains L. Write c_k for the optimum vertical center density in the corresponding lattice or periodic class. The following bound holds either for lattice stacks where L is the full horizontal intersection, or for arbitrary periodic stacks with at most 3^k distinct phases:
\[
\delta\le \frac{c_k}{2^n}\frac{\mu^{(n+k)/2}}{\max\{\mu/9, \mu-D\}^{k/2}}.
\]
For a lattice, the projected lattice has a phase homomorphism and its triple has zero phase, giving projected minimum squared at least μ/9. For a periodic packing, in each fixed phase the vertical minimum is at least μ, and summing at most 3^k classes gives the same density bound. When μ>D each vertical site houses at most one phase and different sites have distance squared at least μ−D, giving the competing estimate in both cases. These estimates do not require equal layer spacing.

The expression increases until μ=9D/8 and beyond that its maximum is at one of the two endpoints 9D/8,D+2: its logarithmic derivative on the second branch changes sign at most once, from negative to positive. Squared ratio of its value at the crossover to its inherited-minimum value is (9D/[8(D+2)])^n*(18/(D+2))^k; the verifier checks it is at most one for k≤12 in 48 and k≤15 in 72. Thus in particular there is no improvement over the established half-period record constructions in 49–53 and 73–76 for these lattices, or in 49–51 and 73–75 for periodic palettes at most 3^k, at **any** μ. The third-period families themselves need not attain the upper bound.

### Arbitrarily large palettes in one vertical dimension

At μ≥D−2/3+10^−4, the same comparison with the established 49D and 73D half-period stacks holds even for arbitrary ternary palettes. The character vectors (exp(2πi x·v/3)) indexed by x on the parent shortest shell have Gram entries N on the diagonal and the coefficient in the table off-diagonal, whenever a phase difference is at its depth cap. Pair depth >D−2/3 forces each attained depth to be the maximum in its own residue because consecutive allowable squared depths in a residue differ by 2/3. The Gram matrix applied to the all-ones vector bounds the number of such phases by M=1+floor(N/a), where a is the smallest magnitude of a negative coefficient in the relevant residues. This gives M=8879 in 48 and 12944 in 72. If pair depths instead exceed D−2/9, only r=0 is possible, giving the sharper M=2687 and 4864 respectively.

If an interval of vertical positions has length sqrt(μ−D+b), with b either 2/3 or 2/9 and μ>D−b, all phase pairs in a half-open interval have depths >D−b and no phase repeats. Its occupancy is at most M, including multiple phases at a height. The squared ratio of its resulting global density bound to the known half-period stack of squared min D+2 and gap sqrt(2) is
\[
2M^2\frac{(\mu/(D+2))^{n+1}}{\mu-D+b}.
\]
Its log derivative changes sign at most once, from negative to positive, so checking endpoints suffices. Exact comparisons cover from D−2/3+10^−4 to D using b=2/3, from D to D+1/12 using b=2/9, and the direct projected distance spacing covers D+1/12 to D+2. Smaller μ for unlimited palettes is not ruled out.

## Fifth-torsion paired theorem

For any v in L define d_i=dist(i v/5,L)^2 for i=1,2. For isotropic v (v²/2≡0 mod5),
\[
\min(d_1,d_2)\le D,\qquad
\min(d_1,d_2)=D\ \Longrightarrow\ \max(d_1,d_2)\le D+4/5.
\]
For any v, d1+d2≤11n/60. If v is anisotropic, this improves to 42/5 in 48D, 64/5 in 72D, and 86/5 in conditional 96D. These are paired statements: alone they do not bound the depth of an individual fifth-torsion shift by D, and do not exclude all fifth-period vertical stacking.

### Scalar coefficients

Put f_v=Σ_x exp(2πi x·v/5)q^(x²/2). For a,b in (1/5)L, the symmetrized theta
\[
G_{a,b}=e^{-\pi i a\cdot b}\sum_{x\in L+a} e^{2\pi i x\cdot b}q^{x^2/2}.
\]
transforms under the weight-n/2 slash action with characteristic row (a,b) acted on by T and S. This follows by completing exponents under T and Poisson Fourier transform under S; (-i)^(n/2)=1. Equivalently T sends (a,b) to (a,a+b) and S sends it to (b,−a). For γ in Γ0(5), the image of G_(0,iv/5) is
\[
e^{-\pi i c d i^2 v^2/25} f_{d i v},
\]
where c,d are the bottom row of γ. Since c≡0 mod5, this scalar is 1 for isotropic v; for the product of i=1,2 it is always 1 because 1²+2²=5. The units modulo 5 permute the two indices up to sign, and f_{−v}=f_v. Thus f_v+f_(2v) in the isotropic case and f_v f_(2v) always are holomorphic scalar forms on Γ0(5). At cusp zero they give the corresponding sums or products of shifted theta series with positive coefficients. Their cusp orders at zero (local variable q^(1/5)) are (5/2) min(d1,d2) and (5/2)(d1+d2), respectively.

Let U=eta(τ)^5/eta(5τ), V=eta(5τ)^5/eta(τ). The eta quotient criterion gives weight two and character (5/d) for both (the two mod 24 sums of exponents are 0 and 24, interchanged); their cusp orders are (0,1) and (1,0) at infinity and zero respectively. For weight 2J with J even, U^j V^(J−j), 0≤j≤J, span the trivial-character space by distinct orders and the index-6 valence upper bound. Here J=6s for the sum and 12s for the product. Their first s+1 coefficients are (2,0,…,0) or (1,0,…,0), and exact normalization makes the coefficient at j=J−s nonzero in both. Therefore cusp zero order is ≤J−s. At equality the leading coefficient is the unique triangular coefficient at j times 5^((j−5(J−j))/2): the factors of U,V under S are −sqrt(5),−5^(-5/2), and J is even. The replay obtains combined counts 60,3060,175200,10563300 at isotropic equality and products 7500,104062500,1623632812500,26685278320312500 at product equality.

### Harmonic and integrality sharpening

For the harmonic theta series of degree h, the same T/Poisson law holds (use the harmonic Gaussian Fourier identity and oriented characteristics). Work on Γ1(5), where the remaining character is a finite fifth root; the possible odd sign is handled by raising the form to its tenth power. The image of Γ1(5) in PSL has index 12. The two Γ0(5) cusps split into two each because the stabilizer T at infinity and ST⁵S⁻¹ at zero lie in the subgroup with bottom-right ≡±1. Thus the first two cusp orders are both at least s+1 for any h>0, while the two width-five cusp orders total at least (5/2)(d_1+d_2). Their sum greater than weight 12s+h forces vanishing. In particular:
* If the isotropic min is D and max≥D+6/5, the total minimum order is at least 12s+5, so h=2,4 vanish. The nearest shell at depth D has exactly the full combined count M from above (the other phase is deeper). Taking x in that shell, each x·y is an integer, since x·y=D−(x−y)^2/2. The spherical fourth moment would give Σ_y (x·y)^4=3 M D^4/[n(n+2)], not an integer in any rank tested. Values are 60/13,4896/5,4730400/37,676051200/49. Allowed isotropic depths have step 2/5, so max≤D+4/5.
* If the anisotropic product cap were attained, orders total at least 13s+2, forcing every degree 1 to s+1 to vanish. For a nearest shell of squared radius q, size M, and a fixed shell point x, all q−x·y are integers. Its (s+1)-st shifted moment must equal M q^(s+1) times the factor shown here: (n+3)/n for s=2; (n+3)(n+5)/(n(n+2)) for s=3; (n+5)(n+7)/(n(n+2)) for s=4. This follows by expanding and using the spherical first through fifth moments. Each q has denominator 25 in the anisotropic case, so the product of the two shell sizes has 5-adic valuation at least 12,12,20 in those ranks. The exact product values have valuations 7,11,14 respectively, contradictions. The sum d1+d2 has step 2/5 because residues of the two exponents sum to an integer in the cusp variable. This yields the refined bounds.

## Replay and references

Run **python3 verify.py** (standard library, exact). References: [Nebe, rank-48 extremal lattices](https://arxiv.org/abs/1312.3780); [Nebe, extremal rank 72](https://arxiv.org/abs/1008.2862); [Ono, eta quotient transformation and cusp criterion](https://uva.theopenscholar.com/files/ken-ono/files/011.pdf); [48D half-period construction and exact density](https://github.com/alexlegeartis/SpherePackings/tree/main/research/p48n_hurwitz_20261122); [existing 73–76 density table](https://github.com/alexlegeartis/SpherePackings). The earlier comparison scope is also reproduced in the [quarter-period packet](../../scoped_density_obstructions_round15_20261123/48-72-quarter-torus/README.md).
