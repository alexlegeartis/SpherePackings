# The quarter-period depth bound for extremal 48-dimensional lattices

**Theorem.** Let L be any even unimodular lattice of rank 48 and squared
minimum 6. Write `d(t) = min_{x in L+t} ||x||²`. Then `d(t) <= 4` whenever
`4t in L`. For a half-period `t=v/2` with `v² = 2 mod 4` the sharper bound
is `d(t) <= 7/2`. A half-period of depth 4 has exactly 2400 nearest vectors,
and their shell is a spherical 2-design. A half-period of depth 7/2 has
exactly 108 nearest vectors. No assertion about arbitrary real shifts is
made. The result requires no full shell enumeration or classification of
the extremal lattices.

The input minimum and examples can be found in [Nebe, *A fourth extremal
even unimodular lattice of dimension 48*](https://arxiv.org/abs/1312.3780)
and the [Nebe–Sloane catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ48a.html).
The Gauss-sum/Poisson identities below follow by Fourier transform of a
translated Gaussian; for ordinary and harmonic theta background see
[Miezaki–Nebe](https://arxiv.org/abs/2009.02842). All short Taylor jets and
the one moment computation replay exactly with `python3 verify.py`.

## 1. Twining formulas and two theta rings

Set `q = exp(pi i tau)` and `Q=q²`. Define
`Psi_(a,b)(tau) = sum_{x in L+a} exp(2 pi i (x,b)) q^{x²}`. With the
weight-24 slash operator, Poisson summation and evenness give

```
Psi_(a,b)|S = exp(2 pi i (a,b)) Psi_(b,-a),
Psi_(a,b)|T = exp(-pi i a²) Psi_(a,b+a).
```

Use the twining `Phi_t=Psi_(0,t)`. It has jet `1+O(Q³)` since L has squared
minimum 6. Its transform by S is the **untwisted** theta of `L+t` with
positive coefficients, so its order at cusp 0 of width w is `w d(t)/2`.
If `mt=v in L`, conjugation by `S T^m S` gives a scalar
`exp(pi i m t²)`. For m=2,4, Γ0(m) is generated projectively by this
parabolic and T; weight24 ignores −I. Thus there is a scalar character on
the whole group. For m=4 this also follows by conjugating Γ0(4) to Γ(2)
by the change `tau -> 2 tau`.

At level 2 set `F=(theta3(tau) theta4(tau))²`, `U=theta3(tau)^4+
theta4(tau)^4`, and `V=F²`. `U,V` are trivial-character forms of weights
2,4, nonzero at infinity; `U` is nonzero at 0, and `V` has order 1 at 0.
`F` has weight 2, sign − under `S T² S`, and order 1/2 at 0. Theta
constants are nonzero on the upper half plane. The usual identities under
S and T verify these statements (e.g., `F|S=-(theta3 theta2)^2` and
`U|S=-(theta3^4+theta2^4)`). The valence total of a nonzero weight-k form
on Γ0(2) is `k/4`. The indicated monomials in U,V below span the residual
weights by their independence and this valence bound; no unproved dimension
assumption is needed.

At level 4 set `A=theta2(2 tau)^4`, `B=theta4(2 tau)^4`,
`C=theta3(2 tau)^4=A+B`. They have weight 2 and trivial character on
Γ0(4), with a simple zero respectively at infinity, 0, and 1/2, and none
in the upper half plane or the other cusps. For 0 this follows by S and
scaling, and for 1/2 by the transformed algorithm given below. Frac-
tional powers used below mean powers of the nonvanishing theta constants
with the branch taking leading value 1 at infinity for B,C. After dividing
by the displayed B,C powers the residual character is trivial: T has
phase 1, and at the other generating parabolic the phase is the matching
fractional cusp order. The cusp orders show the quotient is holomorphic.
The valence total on Γ0(4) is `k/2`; it also follows from its index 6 in
SL2(Z). The forms `A^j B^(k-j)` span residual weight 2k because their
orders 0,...,k at infinity are different and the valence bound gives
dimension at most k+1. Use this argument also for k=0.

## 2. All half-periods, including contact counts

Suppose `t=v/2` is not 0 modulo L. The depths are integers for `v²=0 mod4`
(the singular case) and half-integers for `v²=2 mod4`; all are at least
`6/4`. In the singular case `d>=5` makes `Phi_t` a combination of `V^5 U²`
and `V^6`. Their first three coefficients in Q are `(4,-128,-384)` and
`(1,-96,4512)`. After imposing constant 1 and linear coefficient 0, the
quadratic coefficient is −2400, a contradiction. For the nonsingular
case `d>=9/2` forces `F^9 U³` and `F^11 U`, with coefficients
`(8,0,-6912)` and `(2,-128,3392)`. The same argument gives −864. These
prove the two bounds.

At singular equality divide by `V^4`; the remaining weight-8 monomials
are `U^4,U²V,V²`. The infinity jet forces coefficients
`(75/2048,21/256,11/128)`. After S only the first contributes to the
depth-4 coefficient, giving `(75/2048)*16^4=2400`. At nonsingular depth
7/2 divide by `F V^3`; the residual weight-10 monomials are `U^5,U³V,UV²`,
with coefficients `(27/4096,37/512,27/256)`. Under S the signs on U and F
cancel, giving `(27/4096)*4*16^3=108`. In each computation the exact
forms follow from the valence dimension bound and the invertible 3-jet
system in `verify.py`.

For a degree-2 harmonic polynomial P, insert P into `Phi_t`; the same
Poisson formulas hold in weight 26, up to a constant unit factor for S
which does not affect orders. The squared-depth-4 half-period is singular,
so the character at Γ0(2) is trivial as above (the two S factors cancel).
The harmonic form has orders at least 3 at infinity and 4 at 0. This
exceeds `26/4`, so it is zero. Every nearest shell has zero harmonic
quadratic sums and is antipodal, hence is a 2-design. This argument is
about the equality shell at depth 4, not a covering-radius proof for
all real phases.

## 3. Norm restrictions at the third cusp for a quarter-period

Put `t=r/4`, `h=2t`, and assume h is not zero modulo L. At cusp 1/2 apply
`gamma=S T^(-2) S`, projectively `(1,0;2,1)`. The transformation is
`Phi_t|gamma = exp(-2 pi i t²) Psi_(2t,t)`. Indeed the two steps above
first yield `Psi_(-2t,-t)` times that scalar, which equals `Psi_(2t,t)` by
negation. Because `3t = -t+r`, replacing the twist `-t` by `3t` gives
factor `exp(2 pi i (2t,r))=1`; here `r²` is even. Negation fixes the shift
coset `L+2t`. Thus `Psi_(2t,t)|T=exp(-pi i r²/4) Psi_(2t,t)`.
Only norms congruent to `-r²/4 mod2` can have a nonzero coefficient at
this cusp. At any specified norm its coefficient has absolute value at
most the number of vectors of that norm in `L+h`: the summands have
modulus 1. **Equality of this magnitude with the unweighted count is
not assumed.** A cancellation only helps the upper bound.

Assume `d(t)>4` for contradiction. The fractional part of a squared norm
of `L+t` modulo 1/2 is `r²/16`; the actual residues imply the lower
orders a for cusp 0 shown below. For the third cusp, b is a lower order
using the just-proved congruence and the known half-coset shell.

| r² mod 8 | d(h) | a | b | absolute coefficient bound at order b if used |
|---:|---:|---:|---:|---:|
|0|2|9|1|2|
|0|3 or 4|9|2|not needed|
|4|2 or 4|17/2|5/2|not needed|
|4|3|17/2|3/2|96|
|2|5/2 or 7/2|33/4|7/4|100 or 108|
|6|5/2|35/4|5/4|2|
|6|7/2|35/4|9/4|not needed|

Here is the elementary shell check. For two numerator vectors v,w in the
same class modulo 2L, with `v != ±w`, the two nonzero lattice vectors
`(v+w)/2,(v-w)/2` both have norm at least 6. Their norms sum to
`(v²+w²)/2`. Consequently a norm-8 or norm-10 shell in such a class is
one antipodal pair, and norm 8 and norm 12 cannot coexist. Distinct
norm-12 antipodal pairs are mutually orthogonal, hence there are at most
96 signed vectors of norm 12 in their class. This proves all the counts
except the possible norm-14 shell in a class of depth 5/2 and the
already proved count 108 at depth 7/2.

For the remaining count take its fixed norm-10 numerator v. Every
norm-14 vector w in that class is orthogonal to v by the same sum argument,
and maps bijectively to a norm-6 vector `x=(v+w)/2` with `(x,v)=5`.
Conversely `w=2x-v` verifies all the properties. All norm-6 vectors
have `|(x,v)|<=5`: a larger integer would give a nonzero vector `v±x`
of norm ≤4. The norm-6 shell has 52,416,000 elements (this follows from
the weight-24 theta series with jet `1+O(Q³)`, also checked in the script)
and is an 11-design. Indeed each even harmonic weighted theta of degrees
2 through 10 is divisible by `Delta³` but has weight below 36, hence is
zero; odd degrees vanish by negation. Apply those spherical moments to
the degree-10 polynomial in `(x,v)` which is zero at absolute values
0,...,4 and 1 at absolute value 5:
`prod_{i=0}^4 ((x,v)^2-i^2) / prod_{i=0}^4 (25-i^2)`.
The summed value is exactly 200. Half lie at +5, so the desired count is
100. [Nebe–Sloane's P48n entry](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ48a.html)
also records the shell size.

The missing depth `d(h)=3/2` has a short direct proof. Take a norm-6
vector v with `h=v/2 mod L`, so `t=v/4+w/2` up to L. A closest vector
`x` in the half-coset `L+w/2` has an antipode in the same coset. Hence
`d(t) <= ||x||² + 3/8 - |(x,v)|/2`. If `d(w/2) <= 7/2` this is <4.
Otherwise its depth is 4 and its nearest shell is a 2-design. The mean
square of `(x,v)` on that shell is `4*6/48=1/2`; values are integers or
half-integers, so there exists an absolute value at least 1. The same
inequality then gives depth less than 4. This exhausts the table inputs.

## 4. All quarter cases contradicted by the 3-jet

In every row of the table divide `Phi_t` by `B^a C^b`. As justified above,
the quotient is holomorphic of weight `24-2(a+b)` and trivial character.
If `a+b=11`, it is `B+u A` after enforcing constant coefficient 1. The
coefficient of Q forces u; the Q² coefficient is nonzero in every case:

| a | b | u | Q² coefficient |
|---:|---:|---:|---:|
|9|2|4|−1632|
|17/2|5/2|7/2|−1216|
|35/4|9/4|15/4|−1416|

If `a+b=10`, the quotient is `B²+u AB+v A²`. Its three coefficients are
forced by `1+O(Q³)`. At cusp 1/2, the change `tau -> 2 tau` sends gamma
to `T S T`. The weight-2 transforms are `A->C`, `B->-B`, `C->A`, by
the theta identities. Thus the leading coefficient there has magnitude
`16^b |1-u+v|`, since `A=16Q+O(Q³)` and B,C have constant 1. Any
fractional-power scalar has modulus 1 and does not affect this magnitude.

| a | b | u | v | forced magnitude | actual upper bound |
|---:|---:|---:|---:|---:|---:|
|9|1|5|83/8|102|2|
|17/2|3/2|9/2|33/4|304|96|
|33/4|7/4|17/4|233/32|516|108|
|35/4|5/4|19/4|297/32|177|2|

Each forced coefficient is impossible. The proof excludes *all* greater
depths by vanishing-at-least orders, not just equality at the first
candidate. This proves the theorem.

## 5. Consequence for 49–53D full fibres; a concrete remaining target

Fix the untwisted L as period, and place full translates with phase from
`(1/4)L/L` at arbitrary sites in `R^r`. If the full squared minimum is mu,
with `33/8<=mu<=6`, then each pair of vertical sites is separated by at
least `sqrt(mu-4)` by the theorem (coincidence is impossible). For r=1,2,3
use the exact unrestricted normal packing constants (interval line,
hexagonal, Kepler). For r=4,5 use the exact Hermite constants **only when
the normal sites are a lattice**, as they are in a lattice lift. These
constants are recorded at [Cohn's packing tables](https://cohn.mit.edu/sphere-packing/)
and the 3D result is [Flyspeck](https://doi.org/10.1017/fmp.2017.1).
Relative to the live minimum-6 normal constructions the center-density
ratio is at most

```
(mu/6)^((48+r)/2) (2/(mu-4))^(r/2).
```

The log derivative has at most one critical point, a minimum, so it
suffices to check both endpoints. At mu=6 this is 1. At mu=33/8 the
worst r is 5, and its squared ratio is `(11/16)^53*16^5<1`.
The comparisons are the [P48q 49D deep binary lift](https://github.com/alexlegeartis/SpherePackings/tree/main/research/covering_and_deep_cosets_20261115),
the [50D triangular deep-plane lift](https://github.com/alexlegeartis/SpherePackings/tree/main/research/adaptive_layers_20261115),
and the [51–53D P48n Hurwitz lifts](https://github.com/alexlegeartis/SpherePackings/tree/main/research/p48n_hurwitz_20261122).
The exact reference density is `(3/2)^((48+r)/2)` divided by the covolume
of the respective minimum-2 normal `A1,A2,A3,D4,D5`, namely
`sqrt(2),sqrt(3),2,2,2`. The theorem neither bounds low-radius fibres nor
asserts a global bound for nonlattice normal packings in ranks 4 and 5.

The constructive test left over is precise: a real phase of squared depth
`d>4` immediately beats the live 49D lower bound **by a lattice**. Adjoin
the single generator `(t,sqrt(max(6-d,3/2)))` above L. Neighbouring layers
meet the minimum-6 constraint by definition, and all layers two or more
steps apart are separated vertically by at least sqrt(6), regardless of
their horizontal shift. This is a lattice even when t is irrational. Its
gain factor over the live exact 49D value is `sqrt(2/max(6-d,3/2))`.
Conversely if the squared covering radius of L is at most 4, every
minimum-6 packing of full real translates of L in 49D has normal gaps
at least sqrt(2), so none can beat the live value. Thus the exact
minimum-6 full-fibre question at fixed L is equivalent to its squared
covering radius being greater than 4. The present theorem says that such
a phase cannot have order dividing four. A deeper phase has not been
certified, so no density table is changed.
