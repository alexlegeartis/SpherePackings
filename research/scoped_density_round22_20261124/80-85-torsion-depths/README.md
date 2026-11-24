# Rank-80 torsion depth advances and an exact fixed-period 85D packing obstruction

No new unrestricted packing density; no exact attainment claimed for any new depth upper bound. Let L be **any** even unimodular lattice of rank 80 with minimum squared norm 8, e.g. the [Bachoc–Nebe construction](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf). Write d(v/p) for the **squared** distance to L. The input half-shift caps are 7 in the integral class and 13/2 in the half-integral class; the depth-7 shell is an antipodal 2-design. The input third-shift caps for residues v²/2=0,1,2 mod3 are (20/3,62/9,58/9). These are established in sibling bn80_torus85_exact_20261123; that separate arithmetic proof and formulas are reused unchanged.

## Universal rank-80 fourth-torsion bound

For every v in L,
```
d(v/4) <= 29/4.
```
The scalar method extends the exact 48/72 proof (sibling quarter_torus_exact48_72_20261124, PRIOR_GATES.md and README.md). The details of the extension are as follows. Set H=d(v/2), r=v²/2 mod4. H=0 reduces immediately to a half phase. Otherwise H>=2 and H belongs to r/2+Z, with the above half upper bounds. Pick u=v−2w with u²=4H. Then v/4=w/2+u/4. If z is a nearest point of w/2+L of squared norm D0, antipodality gives
```
d(v/4) <= D0+H/4−max |z·u|/2.
```
Here z·u belongs to Z+c with c=((r−2H) mod4)/4 in {0,1/2}. When D0<7 use D0<=13/2, |z·u|>=c; when D0=7 the design gives the mean squared dot 7u²/80=(7/20)H, so use the least nonnegative s in Z+c with s²≥(7/20)H. Round down to the progression r/8+(1/2)Z. These comparisons are exact in quarter80.py.

For the quartic twist bring Γ0(4) to Γ(2) by z=2τ. Put X=theta2(2τ)^4, Y=theta4(2τ)^4, W=X+Y. Divide by P=theta4(2τ)^r theta3(2τ)^(4−r) if r>0 (P=1 otherwise). This has the twist multiplier, no interior zeros, and weight2. Its powers at the quarter and half cusps are respectively r/16 and α=(4−r)/4 in Q=exp(2πiτ); both zero if r=0. The residual homogeneous form has degree 20−[r>0]. A proposed quarter depth E requires factor Y^i for i=ceil(2(E−r/8)), and H at the half cusp requires factor W^j with j=ceil(H/2−α). These factors are relatively prime. The remaining homogeneous degree is d=20−[r>0]−i−j. Determine whether coefficient prefix (1,0,0,0) at infinity is possible in P Y^i W^j times the full monomial basis of degree d. There are ten surviving rows above depth7, and only three with E>29/4.

The half cusp M=ST^(−2)S⁻¹ is the twist of the shifted half coset v/2+L, so the complex coefficient at Q^(s/2) has modulus at most the actual unsigned shell at squared radius s (by Poisson). On the level-4 theta basis, (X,Y,W) goes to (W,−Y,X) and P becomes theta4^r theta2^(4−r) up to a common unit. Each of the three rows above29/4 has a unique solution; the coefficient modulus at s=2(α+j) is listed by exact expansion.

| Residue r | H | Hypothetical quarter E | s | Required modulus | Unsigned cap |
|---|---:|---:|---:|---:|---:|
| 0 | 4 | 15/2 | 4 | 9,920 | 160 |
| 3 | 7/2 or 9/2 | 59/8 | 9/2 | 17,572 | 12,800 |

At s=4 non-antipodal half-coset vectors are orthogonal (x±y lie in L and have norm²≥8). At s=9/2 they have absolute dot exactly 1/2 by the same bound and half-integral quantization. The corresponding lines have coherence 1/9. Their trace/rank bound is 80(1−1/81)/(1−80/81)=6,400 lines. This excludes all higher depths. The table proves a bound on fourth torsion only.

## Fifth-torsion paired and pointwise bounds

For d1=d(v/5), d2=d(2v/5),
```
d1+d2 <= 72/5  (all residues)
min(d1,d2) <= 34/5 if v²/2=0 mod5
max(d1,d2) <= 202/25; <=8 if v²/2=0 mod5.
```
The fifth-torsion scalar construction extends the earlier sibling third_torus_extremal_20261124 paired proof. Let f_i be the twists by i v/5, i=1,2. On Γ0(5), f1 f2 is a scalar form of weight80, and f1+f2 is scalar of weight40 in the isotropic case. To verify the multiplier, use the symmetrized characteristic G_(a,b), which transforms by row action on (a,b), and track the factor exp(−πi c d i² v²/25). This is 1 in the isotropic case or cancels in the product because c is divisible by5 and 1²+2²=5. The units permute i=1,2 up to sign. At cusp zero the norm counts in the shifted theta expansions are nonnegative. The two eta weight-two forms U=η(τ)^5/η(5τ), V=η(5τ)^5/η(τ) have respective infinity/zero local cusp orders (0,1) and (1,0) and the same quadratic character. In the even-weight exponent cases J=20,40, the monomials U^i V^(J−i) are a full basis by valence and independence. The leading cusp-zero order is bounded by J−3 using the fixed four infinity coefficients. At equality the product of nearest-shell sizes is 56,591,796,875,000. Scalar provisional bounds are 34/5 for the isotropic smaller depth and 74/5 for the sum.

To exclude sum equality, take harmonic twists of degrees h=1,…,4. On Γ1(5) the two infinity cusps have orders at least4 each and the two width-five cusps have total order at least (5/2)(d1+d2)=37. The two split pairs of cusps follow from the stabilizers T and ST^5S⁻¹ in Γ1(5) up to sign; its PSL index is12. Character phases (including odd parity) can be removed by the tenth power. Since 45>40+h, valence forces vanishing, so both nearest shells are spherical 4-designs (no antipodality needed). For any fixed vector x in a nearest shell of squared radius d and size M, d−x·y is integral. Its fourth moment over the shell is
```
M d^4 (83*85)/(80*82).
```
In the anisotropic case both d have reduced denominator25, so each M is divisible by 5^8, whereas their fixed product has v5=14. In the isotropic case both 5d are positive even integers summing to74, hence neither is divisible by41. The fourth moment then forces each M divisible by41; their product is not. Thus 74/5 is impossible, and sum depths lie in (2/5)Z.

For the individual bound let x be nearest to zero in 2v/5+L, so v/5=x/2+w/2 for w∈L. Using antipodal half contacts gives d1≤13/2+d2/4 if d(w/2)<7, or d1≤7+d2/4−sqrt(7 d2/80) by the 2-design otherwise. If d1>202/25, the paired theorem gives d2<158/25; the first branch contradicts directly, the second is convex in d2 with both endpoint values below202/25. Isotropic quantization rounds down further to8.

## Exact 85D packing limits with fixed pure vertical period 2D5

Take a period lattice L×(2D5), and one full horizontal translate per selected vertical site of its third grid (2/3)D5/(2D5). Horizontal phase differences are in (1/3)L/L. A motif may choose any subset of vertical positions and arbitrary independently assigned phase labels; distinct phase labels are not required. Identify vertical coordinates with F3^5. For a difference of Hamming weight w, the squared torus distance vertically is (4/9)(w+3[w odd]): rounding to entries 0,±1 has even parity if w even, and if w odd a ±3 correction on a nonzero coordinate costs exactly3. Thus every weight-two pair with μ>70/9 is incompatible from the universal third cap62/9. In this range multiple phases at one vertical position are also automatically impossible. The distance-two Cayley graph has order243, degree40, and eigenvalues (40,16,1,−5,−2,10) by Fourier character weights. Its independence bound is 27. The 85D incumbent with μ=8 corresponds to 32 points per this same period volume, using half phases. Therefore no third palette on this period improves when μ>70/9.

More generally no palette of up to 81 motifs beats the incumbent at **any** μ. For μ≤70/9 the squared density ratio is at most (81/32)² (35/36)^85 <1. This closes every additive subgroup of rank≤4 with this fixed pure vertical period. For lower μ an arbitrary single-site-per-class palette of ≤243 also cannot improve when μ≤8(32/243)^(2/85), which is strictly above68/9. General repeat occupancy at low μ and other vertical periods/scales are not covered.

A full rank-five additive improvement would therefore require 68/9<μ≤70/9 together with the following stringent conditions. Every Hamming-weight-two phase must have residue1 and attain squared depth62/9, because all other residue caps or the next largest admissible r=1 depth give lower μ. For a linear phase map φ, the quadratic form Q(φ(w)) mod3 must be exactly 2||w||_H mod3: testing e_i±e_j forces zero polar cross terms and equal unit-axis norm coefficients2. Consequently all 90 weight-one and weight-four classes together (10+80) must have r=2 and attain 58/9, since their vertical cost is16/9 and the next possible depth is52/9. These constraints are necessary, not a construction and not an impossibility proof for the rank-five window.

## Replay

Run `python -B quarter80.py`, `python -B fifth80.py`, `python -B palette85.py` (Python 3 and SymPy). All theorem comparisons are exact. Rank80 half and third inputs: bn80_torus85_exact_20261123/verify.py. Quarter and fifth eta/theta transformation derivations are in the named prior sibling packages. Eta quotients use the standard cusp criterion (see [Ono](https://uva.theopenscholar.com/files/ken-ono/files/011.pdf)). No Git changes.
