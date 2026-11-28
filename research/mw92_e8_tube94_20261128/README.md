# An eight-dimensional tube cover improves the 94D MW92 lattice packing

There is a 94-dimensional **lattice** packing of occupied fraction at least **8.08530552528770e−21**, more than **0.0874%** above the maintained [94D paired-cover bound](../mw92_paired_cover_9495_20261126/), whose occupied fraction was at least \(8.078244769088839\cdot10^{-21}\). Both the lattice-only and unrestricted lower bounds therefore improve relative to that named benchmark. As in the predecessor, the horizontal shears are existential; this is not a numerical 94-by-94 Gram matrix or an optimality claim.

The new ingredient is geometric. A sparse covering of a literal scaled \(E_8\) plane covers **whole tubular neighborhoods** of the weak phase. The replacement costs fewer large balls than separate dimers and uses an ordinary periodic equal-ball covering theorem, not an inequality for mixed radii. An elementary \(A_8^*\) version already beats the old bound; Schürmann and Vallentin's published eight-dimensional covering gives the displayed stronger result.

## Parent, exact rank-eight lattice, and benchmark

Use [Shioda's characteristic-47 construction](https://doi.org/10.2307/2374791) in the same normalization as the [audited MW92 parent](../mw92_extensions_20261122/): the even integral lattice \(L=\mathrm{MW}_{92}\) has squared minimum 16 and covolume \(47^7\). The elliptic model is
\[
y^2=x^3+1+t^{48},\qquad k=\mathbb F_{47}[\omega],\quad \omega^2+\omega+1=0.
\]
Multiplication of \(x\) by \(\omega\) yields an integral isometry \(J\) with \(I+J+J^2=0\). All sections below are literal over \(k(t)\); no claim about descent of all geometric sections is needed.

The [nine-site parent certificate](../mw92_nine_cluster95_20261128/) checks eight literal sections \(v_1,\ldots,v_8\) in the variable \(u=t^8\), with Gram \(8(I_8+\mathbf1\mathbf1^{\mathsf T})\). Thus they generate a copy of \(\sqrt8 A_8\), identifying \(v_i/\sqrt8\) with \(e_i-e_9\). In the same model add the literal section with polynomial coefficient arrays
```
x: [864, 1029, 2]
y: [2161, 1920, 245, 44]
```
in ascending powers of \(u\), with field encoding \(a+47b=a+b\omega\). This is the base-inversion transform of the first archived section in the small model, but the checker uses its coefficients directly. It verifies the curve equation and pairings
\[
 \|w\|^2=16,\qquad
 (\langle w,v_i\rangle)_{i=1}^8=(0,0,0,0,8,0,8,8).
\]
These imply \(w=v_5+v_7+v_8-(v_1+\cdots+v_8)/3\): subtract the right hand side and compute its squared norm as zero. Therefore adjoining \(w\) is precisely the index-three glue extension of \(A_8\) by the fundamental weight of a three-element subset. The resulting rank-eight sublattice \(S\subset L\) is isometric to \(\sqrt8 E_8\), with covolume \(8^4\) in its span \(W\). Only inclusion in \(L\) is used. To check heights, the companion script expands \(u=t^8\) in \(k[t]\) and uses the published minimum and the [original model's height computation](../mw92_cusp_d4_20261122/): the 48 bad fibers are irreducible type II, infinity is smooth, degree bounds make the sections disjoint from the zero section, and GCD intersection counts (also at infinity) give their pairings.

The exact previous squared vertical covolume for 94D is \(850355653/45000000\) in the [earlier arithmetic ledger](../mw92_paired_cover_9495_20261126/result.json). Unlike BCS-style derived products, this is an exact comparison to the maintained 94D lower bound.

## A periodic equal-ball replacement for the weak phase

Write \(d(x)=\operatorname{dist}(x,L)\) and use
\[
 A=11.00098526,\quad s=16-A,\quad p=14/3,\quad D=10.79109,\quad r^2=D-2p .
\]
The established [Coxeter–Few–Rogers input](../mw92_cfr_taylor13_20261123/) says that an **ordinary equal-ball** cover of \(\mathbb R^{92}\) of squared radius \(D\) with at most \(N\) centers per \(L\)-cell is impossible if \(N D^{46}<A^{46}\).

Choose inside \(W\) a full rank lattice \(C\) commensurable with \(S\), of covering radius \(\rho\) and covolume such that
\[
 \rho^2\le 2p,\qquad 
 \beta=\operatorname{covol}(S)/\operatorname{covol}(C)<1703/4000.
\]
Here is the precise source and rationale. By [Schürmann and Vallentin, Theorem 1.2 and Appendix A](https://arxiv.org/abs/math/0405441), there exists a rank-eight lattice with covering density below 3.2013, in the convention \(\Theta=V_8\rho^8/\operatorname{covol}(C)\), \(V_8=\pi^4/24\). Scale its covering radius to \(\sqrt{2p}\) initially. Then its ratio to \(S\)'s covolume meets the strict target because
\[
 \frac{3.2013\cdot8^4}{(\pi^4/24)(2p)^4}<\frac{1703}{4000}.
\]
This last comparison is checked with \(\pi>3.14159\) exactly. Now shrink very slightly (so that the radius has strict slack), then approximate that lattice basis by rational linear combinations of a basis of \(S\). Determinant is continuous; more explicitly, if the basis change from the shrunk lattice to the approximation is a real linear map \(M\) tending to identity, the covering radius of the latter is at most \(\|M\|\) times the former. Hence both strict inequalities persist, and the rational coordinates imply commensurability. This uses only the published existence of the eight-dimensional covering, not an approximation as if it were an exact basis.

In any affine translate \(x_0+W\), replace \(r\)-balls centered at \(x_0+S\) by \(\sqrt D\)-balls centered at \(x_0+C\). Every point of the old balls is within perpendicular distance \(r\) of the plane, and its projection has a covering point of \(C\) at squared distance at most \(2p\). Pythagoras gives squared distance at most \(r^2+2p=D\). Covering the whole tube proves every needed point, even if no individual large ball contains a whole small ball.

The replacement is periodic with the stated count. Set \(S_{\mathrm{sat}}=L\cap W\). It is a primitive rank-eight sublattice, with a free complement \(L_1\) in \(L\); put \(H=S\cap C\). For each of the finite classes \(S_{\mathrm{sat}}/S\), represented by \(x_j\), replace the weak centers \(L_1+x_j-u+S\) by \(L_1+x_j-u+C\). The combined centers are periodic under \(L_1+H\subset L\), which is full rank. The weak old and new counts per such cell are \([S_{\mathrm{sat}}:S][S:H]\) and at most \([S_{\mathrm{sat}}:S][C:H]\), respectively. Their quotient is at most \(\beta\); any coincident new centers only improve the count. No saturation assertion for \(S\) itself is needed.

## Exhaustive two-shear construction

The original one-ball covering threshold yields \(a\) with \(d(a)^2>A\). Set \(T=I-J\) and \(t_0=16-3s\).

1. If \(d(Ta)^2\ge t_0\), keep the previous generic pair \(a,-Ja\) and its triangular vertical Gram of scale \(s\). Norm-one and norm-three Eisenstein vectors are protected by the two depth inequalities above; all remaining nonzero Eisenstein norms are at least four and \(4s>16\). Its exact squared vertical covolume is \(3s^2/4=187426112780529507/10000000000000000\).
2. Otherwise, the [established order-three rounding argument](../mw92_paired_cover_9495_20261126/) supplies \(u\) of order three modulo \(L\), fixed by \(J\), with \(d(u)^2\ge34/3\). Briefly, subtract \(T^{-1}\) of a shortest lift of \(Ta\) from \(a\), obtaining a fixed class \(u\) at squared displacement less than \(t_0/3\). The equal-length relative vectors from \(u\) to an affine-rotation orbit of a nearest lattice point sum to zero; averaging them yields \(d(u)^2>A-t_0/3=32/3\). Since \(T^{\mathsf T}T=3I\), the values of \(3d(u)^2\) are even integers, giving the stated quantized lower bound.

In the second case find a horizontal shear \(b\) with
\[
 d(b)^2>D,\quad d(b-u)^2>D,\quad d(b+u)^2>D-2p .
\]
Its nonexistence would give a cover by two strong \(\sqrt D\)-ball cosets at \(L,L+u\) and one weak \(r\)-ball coset at \(L-u\). Replace the weak coset by the tubes above. The resulting ordinary equal-ball cover has at most \(2+1703/4000\) balls per \(L\)-cell, contradicting the exact arithmetic inequality
\[
 (2+1703/4000)D^{46}<A^{46}.
\]
Hence \(b\) exists. The same earlier packing argument now uses \(q=16-D\), the vertical Gram \(\begin{pmatrix}p&p/2\\p/2&q\end{pmatrix}\), and horizontal shear \(i u+j b\). The quadratic form is \(p(i+j/2)^2+(q-p/4)j^2\); \(|j|\ge2\) is safe since \(q-p/4>4\). For \(j=0\), \(|i|=1\) is safe by \(d(u)^2\ge34/3\) and the other \(i\ne0\) are long. For \(j=1\) the pairs \(i=0,-1\) need the two strong inequalities; \(i=1,-2\) use the same weak class modulo \(3u\in L\) with vertical norm \(q+2p\); any remaining \(i\) has vertical norm at least \(q+6p>16\). Negation treats \(j=-1\). Each branch is a single lattice, of attained squared minimum 16 from the embedded parent. Its covolume is \(47^7\) times its vertical covolume.

## Numeric certificate and elementary backup

The worst squared vertical covolume is the fixed branch:
\[
 D_{\mathrm{tail}}=p(q-p/4)=\frac{8488711}{450000}.
\]
An exact occupied-density lower bound is
\[
\frac{2^{94}\pi^{47}}{47!\,47^7\sqrt{8488711/450000}}.
\]
The gain over the pinned preceding construction is \(\sqrt{(850355653/45000000)/(8488711/450000)}-1>0.000874\). The displayed decimal is rounded down by an exact comparison to the prior certified floor.

For a proof that does not depend on the external eight-dimensional covering theorem, take \(C=\sqrt8(64/51)A_8^*\) parallel to the identified \(A_8\) in \(S\). Its squared covering radius is at most \(8(64/51)^2(20/27)<28/3\), and its exact ratio \(3(51/64)^8<1/2\). For the radius bound, let the coordinates of a Voronoi vector of \(A_8^*\) sum to zero and sort them decreasingly. Squared distance to each projected indicator of a size-\(j\) subset implies that the sum of the top \(j\) coordinates is at most \(j(9-j)/18\). Thus the vector is majorized by \((4,3,2,1,0,-1,-2,-3,-4)/9\), whose squared length is \(20/27\). The determinant of \(A_8\) is 9, so its dual covolume is \(1/3\). Using \(D=10.78516\) makes the same exact strict covering contradiction, yielding squared tail covolume \(2125291/112500\), still smaller than the old bound though not as good as the main result.

Run `node verify.js` from this folder. It verifies the literal \(E_8\) glue, all covering ratios, both exact determinants and the conservative decimal floor. The parent minimum, ordinary covering threshold and better eight-dimensional covering are separately cited published inputs.
