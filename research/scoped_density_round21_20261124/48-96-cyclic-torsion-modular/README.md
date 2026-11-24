# Cyclic torus depth obstructions from the modular curve

Let \(L\) be extremal even unimodular of rank \(n=24s\), \(s=1,2,3,4\), with squared minimum \(2s+2\); use \(D=2s\). Rank 96 is conditional on such a parent existing. For \(v\in L\) and exact order \(N\) of \(v/N\bmod L\), call the line isotropic if \(v^2/2\equiv0\bmod N\). All results concern the simultaneous depths of cyclic phases, not an upper bound on individual shifts or on unrestricted sphere packing.

**Results.** In any isotropic cyclic line with positive-genus \(X_0(N)\), at least one nonzero phase has squared depth strictly below \(D\). For prime \(p\) this gives at most \(D-2/p\) when \(p\notin\{2,3,5,7,13\}\). In the five genus-zero prime levels, some nonzero isotropic phase has squared depth at most \(D\); when all reach \(D\), their entire shifted theta series coincide. Moduli 7 and 13 also have new harmonic-arithmetic bounds on the summed depths of anisotropic lines. A separate modulus-six result is \(3d_6+\tfrac32d_3+d_2\le11s\), with \(d_j=\operatorname{dist}(v/j,L)^2\). No density gain found.

## The residue obstruction at positive genus

Write \(w=12s\) and
$$
 f_i(\tau)=\sum_{x\in L} e^{2\pi\mathrm i\,i x\cdot v/N} q^{x^2/2},
 \qquad q=e^{2\pi i\tau}.
$$
For characteristics \(a,b\), use the symmetrized theta \(G_{a,b}=e^{-\pi i a\cdot b}\sum_{x\in L+a}e^{2\pi i x\cdot b}q^{x^2/2}\). The usual weight-\(w\) slash action of \(T,S\) acts on the row \((a,b)\) by the corresponding matrices: complete the square under \(T\) and use Poisson summation under \(S\), where \((-i)^w=1\). Thus for \(\gamma\in\Gamma_0(N)\) of bottom row \((c,d)\),
$$
 f_i|_w\gamma=\exp(-\pi i c d\,i^2v^2/N^2)f_{di}.
$$
The phase is 1 if the line is isotropic. The sum \(F=\sum_i f_i\) over the \(r\) unit classes modulo sign is a scalar form on \(\Gamma_0(N)\), with \(F=r+O(q^{s+1})\). If every nonzero shift has squared depth at least \(D\), then at each cusp other than infinity, the first characteristic \(c i v/N\) is nonzero in the torus, so every summand vanishes to order at least \(s\) times the cusp width. (Cusp infinity is exactly the orbit of denominators divisible by \(N\).) The form may have signed or complex coefficients there, which cannot reduce this order. Thus \(H=F/\Delta^s\) is a meromorphic function on \(X_0(N)\) with only one pole, at infinity of order \(s\).

If \(X_0(N)\) has positive genus there is a normalized weight-two newform on some divisor level whose pullback \(g=\sum_{j\ge1}a_j q^j\) on \(\Gamma_0(N)\) has \(a_1=1\). This is the newform decomposition; fix any complex embedding. The differential \(g\,dq/q\) is holomorphic, including at cusps. The meromorphic differential \(H g\,dq/q\) has only a possible pole at infinity, where its residue is
$$
 r\,[q^s]\left(g(q)\prod_{k\ge1}(1-q^k)^{-24s}\right).
$$
The coefficient of \(a_1\) is respectively \(1,48,2700,161408\), while the sum of absolute bounds from \(j\ge2\) is at most \(0,3,220,14646\). Here \(|a_j|\le d(j)\sqrt{j}\) (at bad primes the elementary newform local factors also obey this), and \(j\le4\). The residue cannot vanish, a contradiction. If \(N=p\) is prime, all nonzero phases are represented by the units, and their depths, under isotropy, lie in \((2/p)\mathbb Z\), giving the discrete refinement. The genus-zero prime list is \(2,3,5,7,13\) by the standard genus formula.

**Gonality refinement.** For prime \(p\) there are only two cusps. Let \(\gamma\) be the complex gonality of \(X_0(p)\). At zero \(F\) is a sum of positive shifted theta series, hence \(H\) has pole order \(\max(0,ps-pd_{\min}/2)\), and at infinity it has pole order \(s\). Consequently if \(\gamma>s\),
$$
 d_{\min}\le D-\frac{2(\gamma-s)}p.
$$
Together with [Abramovich's complex gonality bound](https://arxiv.org/abs/alg-geom/9609012), \(\gamma\ge(7/800)(p+1)\), this gives an explicit uniform positive deficit for sufficiently large primes. For composite levels, the same argument bounds the degree by \(s+\sum_{C\ne\infty}\max(0,w_C(s-\tfrac12 d_C))\), where \(d_C\) is the minimum squared depth of the characteristic shifts at cusp \(C\), and \(w_C\) is its width. This last inequality is valid even if cancellation raises an actual cusp order.

## Prime levels with an eta basis

For \(p\in\{3,5,7,13\}\), set \(t=24/(p-1)\), \(J=24s/(p-1)\). In the scalar trivial-character space of weight \(12s\), a basis is
$$
 E_i=\eta(\tau)^{pJ-ti}\eta(p\tau)^{-J+ti},
 \quad 0\le i\le (p+1)s.
$$
Both eta congruences hold: the orders at infinity and zero are \(i,(p+1)s-i\), and the character is trivial because \(J,t\) are even. These expressions are holomorphic; their distinct orders prove independence, and the valence upper bound gives spanning. For the higher weights below, multiply \(J\) by the relevant positive integer. The first \(s+1\) coefficients at infinity fix the basis amplitudes \(i=0,\ldots,s\) triangularly. The verifier computes nonzero amplitudes at \(i=s\). Hence for the isotropic sum, cusp-zero order is at most \(ps\), proving \(\min d_i\le D\). Level 2 is as in the prior half-period certificate.

**Isotropic equality rigidity for all levels.** If all nonzero phases have squared depth at least \(D\) then \(V=\sum_{i<j}(f_i-f_j)^2\), over units modulo sign, transforms as a scalar form of weight \(24s\). It has order at least \(2s+2\) at infinity, and at least \(2s\) times the width at every other cusp. Their sum exceeds the valence bound, so \(V=0\). On imaginary \(\tau\) the \(f_i\) are real, forcing all the \(f_i\) identical. For the indicated prime levels, their shifted theta series are thus identical at equality. In particular this strengthens the previously separate fifth-period cap: if both depths reach \(D\), **both equal \(D\)**. The nearest count **per unit phase** is:

| n | p=3 | p=5 | p=7 | p=13 |
|--:|--:|--:|--:|--:|
|24|36|30|28|26|
|48|1836|1530|1428|1326|
|72|105084|87600|81760|75920|
|96|6334524|5281650|4929540|4577430|

Indeed the leading coefficient at cusp zero of \(E_s\) is 1: the eta transformations give phase \((-i)^{12s}=1\) and multiplier \(p^{(J-ts)/2}=1\). Divide the sum's coefficient by the number of phases.

## Anisotropic odd-prime product and proper block bounds

For \(p=5,7,13\), the product of all \(r=(p-1)/2\) twists is scalar even for anisotropic phases because \(\sum_{i=1}^{r}i^2=r(r+1)p/6\equiv0\bmod p\). With \(J\) replaced by \(rJ\), the basis above yields raw width-normalized cusp-zero order at most \(U=r(p+1)s-s\). At a raw equality the product of nearest counts is its triangular coefficient \(c_s\) times \(p^{(rJ-ts)/2}\), replayed exactly.

For each harmonic-weighted twist of degree \(h>0\) use \(\Gamma_1(p)\): the same oriented characteristic action follows from Fourier transforming harmonic Gaussians; the remaining scalar is a finite character and may be removed by raising to its \(2p\)th power. The image in \({\rm PSL}_2\) has index \((p+1)r\). Each of the two cusps of \(\Gamma_0(p)\) splits into \(r\) cusps with unchanged width, because its stabilizer maps trivially to \((\mathbb Z/p)^\times/\{\pm1\}\). Thus the sum of cusp orders is at least \(r(s+1)+U\) at raw equality, so all harmonic degrees satisfying
$$
 r(s+1)+U > \frac{(p+1)r}{12}(12s+h)
$$
vanish individually. Each nearest shell is then a spherical \(h\)-design.

Taking any \(x\) in a shell of radius squared \(d\), size \(M\), the values \(d-x\cdot y\) are integers for other shell members \(y\) (including \(y=x\)). The \(h\)th moment is \(M d^h \gamma_h\) where, for \(h=1,2,3,4,5\),
$$
 \gamma_h =1,\ \frac{n+1}{n},\ \frac{n+3}{n},\
 \frac{(n+3)(n+5)}{n(n+2)},\
 \frac{(n+5)(n+7)}{n(n+2)}
$$
respectively. These follow immediately by expanding and taking the normalized spherical first through fifth moments. In the anisotropic case every \(d_i\) has reduced denominator \(p^2\). Therefore \(v_p(M_i)\ge2h-v_p(\gamma_h)\) for every shell. The exact product valuation fails the resulting lower bound in every mode checked except \(p=5,n=24\), excluding raw equality. Because the sum has step \(2/p\), the caps on the sum of squared depths are:

| n | p=5 | p=7 | p=13 |
|--:|--:|--:|--:|
|24|22/5|44/7|164/13|
|48|42/5|90/7|330/13|
|72|64/5|136/7|496/13|
|96|86/5|26|662/13|

For \(p=13\), the unit group modulo sign has order 6; partition into cosets of its unique subgroups of sizes \(b=2\) or 3. The sum of squares of each block is zero modulo 13 (a sum over a nontrivial root-of-unity subgroup), so the block products are simply permuted. Their scalar sum and eta basis with \(bJ\) give some block width order at most \(A=b(p+1)s-s\). If all blocks reach \(A\), the variance of the block products is forced zero by the same two-cusp order argument as before. Thus all block products are identical, the total width order is \((r/b)A\), and each block's product of counts is the sum's exact count divided by \(r/b\). Apply the harmonic argument with that new total width and compare the predicted product valuation **within a block**, \(b(2h-v_p(\gamma_h))\). The verifier excludes equality throughout, so some block has summed squared depth at most \(2(A-1)/13\), both for the pair and triple partitions. In 48D these are 106/13 and 162/13; in 72D they are 160/13 and 244/13. In 24D the pair cap equals 4, so even anisotropically some order-13 phase has squared depth at most 2, with anisotropy forcing strict inequality.

## Mixed modulus six, with signed cusps

Let \(d_j=\operatorname{dist}(v/j,L)^2\), \(j=2,3,6\), and \(r=v^2/2\bmod6\). The single twist \(f_1\) transforms by a finite scalar character on \(\Gamma_0(6)\), since the units are \(\pm1\). Its four cusp orders (0, 1/2, 1/3, infinity), in respective width-normalized variables, are at least \(3d_6,\tfrac32d_3,d_2,0\). The latter two finite cusps may have signed coefficients, but valuation lower bounds suffice. Their order residues modulo 1 are \(r/6,r/3,r/2,0\): the row's first shift is \(cv/6\) when the cusp denominator is \(c=1,2,3\), and the width is \(6/c\).

The cusp-order matrix of \(\eta(\delta\tau)\), \(\delta=1,2,3,6\), is
$$
 \frac1{24}\begin{pmatrix}
 6&3&2&1\\3&6&1&2\\2&1&6&3\\1&2&3&6
 \end{pmatrix}.
$$
Its inverse is integral, and sends all six residue tuples (taking fractional parts) to integral eta exponents. For any target vector of nonnegative cusp orders with these residues and sum \(12s\), the resulting eta quotient is holomorphic and has the same weight and character as the twist. To see character equality without using a trivial-multiplier eta criterion, each eta term has a finite character on this group after the integer total weight is taken; check local multipliers on the positive cusp stabilizers using the prescribed residues. Their quotient character is then trivial since the image of \(\Gamma_0(6)\) is torsion-free, genus zero, with four cusps, hence parabolic generators (and the two original forms agree at \(-I\)). Torsion freeness follows from \(x^2+1=0\) mod 3 and \(x^2\pm x+1=0\) mod 2 having no roots; genus zero follows from index 12 and the four cusp widths.

Fix \(A=\tfrac32d_3\), \(B=d_2\) and put \(W=12s-A-B\). If cusp-zero order exceeds \(W-s\), the form must be in the span of eta quotients with target orders \((W-i,A,B,i)\), \(i=0,\ldots,s-1\). This follows directly from valence and the distinct initial orders at infinity: after canceling their first \(s\) coefficients, any nonzero remainder would have total order exceeding \(12s\). The verifier normalizes the first \(s+1\) coefficients starting from \(1,0,\ldots,0\); the coefficient at \(i=s\) is nonzero for every residue tuple in the allowable half- and third-depth ranges (50,171,364,629 cases). Therefore \(3d_6+A+B\le 11s\). In particular, if \(d_2=d_3=D\) then \(d_6\le D\). If all three equal \(D\), the leading cusp-zero term is the \(i=s\) eta quotient \(\eta(\tau)^{24s}\) times its triangular amplitude, giving exact sixth-period nearest counts \(29,1477,84539,5096373\). Attainment is not asserted.

## Reproducibility

Run `python3 verify.py` (pure integer and rational arithmetic). See [Nebe on extremal 48D lattices](https://arxiv.org/abs/1312.3780), [Nebe on the 72D lattice](https://arxiv.org/abs/1008.2862), [Ono on eta quotient cusp orders and characters](https://uva.theopenscholar.com/files/ken-ono/files/011.pdf), [Deligne's proof of the Weil conjectures](https://doi.org/10.1007/BF02684373), and [Abramovich on gonality](https://arxiv.org/abs/alg-geom/9609012). The separate earlier half and third bounds are in the [fourth-period comparison](../../scoped_density_obstructions_round15_20261123/48-72-quarter-torus/HALF_BOUND.md) and [third/fifth packet](../../scoped_density_round20_20261124/48-96-third-fifth-tori/README.md).

### Finite P48n constructive scout (no gain)

Separately checked 3,328 denominator-5/6/7/13 shifted phases of the existing P48n parent (832 per denominator), using fixed-seed random vectors, several eigenspace projectors of an order-six parent isometry, and integer roundings of known half-period colors. A heuristic CVP search located 3,238 literal blockers at squared distance at most 4; the remaining 90 samples reduce to half periods and obey the established half theorem. The exact replay checks the catalogue-basis transform and every literal distance using integer Gram arithmetic. Thus no sampled phase can improve the ordinary single-shift 49D stack via depth greater than 4. This is a finite sample, not a claim about all phases at any denominator. All non-half sampled phases have **strict** literal blockers below 4; boundary cases in the first floating scan were continued and exactly replayed. Replay with `python3 verify_scout.py`; `p48_scout.json.gz` contains the parent Gram, its unimodular replacement, each sampled phase and its short blocker.
