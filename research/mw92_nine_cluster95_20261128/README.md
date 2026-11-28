# A nine-site MW92 simplex improves the 95D lattice packing

A 95-dimensional **lattice** packing has occupied density at least **2.133434801047687e-21**, a gain strictly greater than **0.193046%** over the prior **2.129324206844210e-21**. The construction keeps the two 94D branches, but changes how the weak balls are grouped when adding the third tail generator. This is a lattice existence proof; it certifies but does not compute the real horizontal shears. No 94D bound changes.

## Literal nine-site parent cluster

Use [Shioda's characteristic-47 construction (1991)](https://doi.org/10.2307/2374791) in the normalization of [Schütt–Shioda, Theorem 13.35 and §13.4.2](https://link.springer.com/book/10.1007/978-981-32-9301-4): \(L=\mathrm{MW}_{92}\) is an even integral lattice with minimum squared norm 16 and covolume \(47^7\). In the elliptic model \(y^2=x^3+1+t^{48}\) over \(\mathbb F_{47}[\omega]\) with \(\omega^2+\omega+1=0\), \(J:(x,y)\mapsto(\omega x,y)\) is an orthogonal integral automorphism and \(I+J+J^2=0\).

Set \(u=t^8\), and let \(v_1,\ldots,v_8\) be the eight sections whose x/y coefficient arrays below are in ascending powers of \(u\). A field element \(a+b\omega\) is written as \(a+47b\), \(0\le a,b<47\).

```
[2,1029,864]       [44,245,1920,2161]
[864,1995,94]      [2161,1880,277,3]
[1363,1995,2160]   [2161,336,245,3]
[2160,1029,29]     [44,277,329,2161]
[29,1995,2]        [2161,40,1781,3]
[94,1029,1363]     [44,1781,7,2161]
[0,987]            [1,0,0,1]
[0,21]             [46,0,0,46]
```

These literal points were identified by the base change \(u=t^8\) from \(y^2=x^3+1+u^6\). The proof uses the arrays directly rather than relying on a search. The checker expands the coefficients, verifies all eight section equations and computes the entire exact height Gram matrix \(H=8(I_8+\mathbf1\mathbf1^\top)\).

For clarity, the [previous height audit](../mw92_cusp_d4_20261122/README.md) applies: the finite roots of \(1+t^{48}\) are simple, so the only bad fibres are irreducible type-II cusps. The chart \(s=1/t,\ x'=s^{16}x,\ y'=s^{24}y\) is smooth at \(s=0\). All coordinates have degree at most 16 in \(x\) and 24 in \(y\), so each section is disjoint from the zero section, including in the infinity chart; its height is 16. The pairing of distinct sections is 8 minus their intersection number. A section cannot pass through a finite singular cusp, since then order \(y^2-x^3\) would be at least two while \(1+t^{48}\) has order one. At smooth fibre points the common vanishing orders of the coordinate differences give the intersection multiplicity. Thus the intersection number is \(\deg\gcd(\Delta x,\Delta y)+\min(16-\deg\Delta x,\ 24-\deg\Delta y)\), taking the order of a zero polynomial to be infinity. This is exactly the calculation in the checker; no approximate inner products are used.

The nine points \(0,v_1,\ldots,v_8\) are therefore a regular simplex of pairwise squared distances 16. A face with \(k+1\) vertices has centroid squared radius \(8k/(k+1)\): translate any vertex to zero and use the \(k\times k\) Gram \(8(I+\mathbf1\mathbf1^\top)\).

### Periodic four- and nine-site partitions

**Four-site face:** \(v_1,v_2,v_3\) are linearly independent modulo \(2L\). If a nonempty subset summed to zero modulo \(2L\), we could change signs to obtain a nonzero vector in \(2L\) with norm 16 (support 1 or 2, use opposite signs for 2) or 32 (support 3, use two plus signs and one minus sign). All are below the minimum squared norm 64 in \(2L\). Send the three independent images to the three distinct nonzero elements of \(\mathbb F_2^2\), and extend the linear map to \(L/2L\). The four vertices have each of the four possible values exactly once, so their translates by the index-four kernel partition \(L\). Their centroid has squared covering radius 6 for the vertices.

**Nine-site simplex:** A coefficient relation among the \(v_i\bmod3L\) gives a coefficient vector in the kernel of \(H\bmod3\) by pairing with the \(v_i\). Since \(8\) is invertible, that kernel consists exactly of constants: \((I+\mathbf1\mathbf1^\top)x=0\) forces \(x_i=-\sum_jx_j\), and constant vectors satisfy the equation because \(9=0\bmod3\). Assign all eight distinct nonzero elements of \(\mathbb F_3^2\) to the eight \(v_i\). Their sum is zero, so every possible relation is respected, whether or not the full constant vector is an actual relation. Extend this well-defined linear map from the span to all of \(L/3L\). The nine vertices give distinct values, so the translates by the index-nine kernel partition \(L\). The centroid squared radius is \(64/9\).

In either case, a coset of weak balls of radius \(r\) is covered by one coset of centroid balls of radius \(\sqrt{8k/(k+1)}+r\), one new ball for each group of \(k+1\) original centres. This constructs a periodic ordinary cover rather than applying an equal-ball theorem to unequal balls.

## How the fixed-order and generic lattices are obtained

We use the same [audited equal-ball threshold](../mw92_cfr_taylor13_20261123/) \(A=11.00098526\) and [ordinary equal-ball Coxeter–Few–Rogers theorem](https://doi.org/10.1112/S0025579300002059) as the [paired-cover 94D/95D proof](../mw92_paired_cover_9495_20261126/). Our unchanged 94D setup follows the [five-site predecessor](../mw92_five_cluster95_20261128/). They assert that a periodic equal-ball covering of \(\mathbb R^{92}\) cannot have density as low as \(\mathrm{vol}(B_{92}(\sqrt A))/47^7\). Use \(d(x)=\operatorname{dist}(x,L)\). We reproduce the exhaustive dichotomy and indicate its inherited exact proofs.

Choose \(a\) with \(d(a)^2>A\), which exists because otherwise one such equal ball per \(L\)-cell would give a covering. Put \(s=16-A\), \(T=I-J\), \(t_0=16-3s\).

* If \(d(Ta)^2\ge t_0\), the generic first two shears are \(a,-Ja\), with vertical norm \(s(i^2+ij+j^2)\). At Eisenstein norm 1 the horizontal distance is \(d(a)\), at norm 3 it is \(d(Ta)\), and otherwise a nonzero Eisenstein norm is \(\ge4\) with \(4s>16\). This works at the exact shear scale \(s\), without the prior deliberate rounding to \(5\).
* Otherwise take \(h=T^{-1}\) of a shortest representative of \(Ta\bmod L\), put \(u_0=a-h\), and note \(Tu_0\in L\), \(\|h\|^2<t_0/3\), using \(T^\top T=3I\). For the three equally near orbit points of a closest lattice point under the affine map \(v\mapsto Jv+Tu_0\), the displacements from \(u_0\) sum to zero. Averaging the squared distances after adding \(h\) gives \(d(u_0)^2>A-t_0/3=32/3\). Also \(3d(u_0)^2\) is an even integer, since \(T(u_0-v)\in L\) for \(v\in L\). Thus \(d(u_0)^2\ge34/3\), and \(u_0\bmod L\) has exact order three by \((I-J^2)T=3I\). Let \(p=14/3\), \(B=10.78402070\), \(q=16-B\). The inherited two-strong-one-weak covering estimate finds \(b\) with \(d(b)^2,d(b-u_0)^2>B\), \(d(b+u_0)^2>B-2p\): otherwise cover by those three kinds of balls, pair the one weak phase using a primitive minimum vector of norm 16 to make dimers, and use \(2+\sqrt{B-2p}<\sqrt B\), \((5/2)B^{46}<A^{46}\). Its vertical Gram is \(\begin{pmatrix}p&p/2\\p/2&q\end{pmatrix}\) and first two shears are \(u_0,b\).

For either case, adjoining an independent third real vertical generator with horizontal shear \(c\) yields the full lattice generated by \(L\times\{0\}\) and the three vertically projected basis vectors with the indicated shears. The exact vertical Grams are:
\[
G_g=\begin{pmatrix}s&s/2&-s/2\\ s/2&s&-s/2\\-s/2&-s/2&16-D_g\end{pmatrix},
\quad D_g=10.5295;
\qquad
G_f=\begin{pmatrix}p&p/2&p/2\\p/2&q&q/2\\p/2&q/2&16-D_f\end{pmatrix},
\quad D_f=10.54714 .
\]
If \(k>0\) and the vertical norm \(Q_G(i,j,k)<16\), only \(k=1\) occurs; the replay verifies this using \(i_l^2\le Q_G(G^{-1})_{ll}\) and checked inverse bounds. Negation covers \(k<0\), and \(k=0\) is the independent 94D argument (also checked). With \(k=1\), every forbidden phase for \(c\) is centred at minus the corresponding two-shear combination, with squared radius \(16-Q_G(i,j,1)\). Take the maximum when labels collide; for fixed-order \(i\) is only needed modulo 3.

* **Generic:** six labels require at most \(D_g\) and six weak labels require exactly \(D_g-2s\). Four-site faces cover the weak balls since \(\sqrt6+\sqrt{D_g-2s}<\sqrt{D_g}\), giving at most \(6+6/4=15/2\) equal-radius balls per \(L\)-cell.
* **Fixed-order:** six labels require at most \(D_f\); label \((1,0)\) requires \(D_f-2p\); labels \((0,1),(1,1),(0,-2),(1,-2)\) require \(D_f-2q\). Use single balls, a minimum-vector dimer, and the nine-site simplices, respectively. The dimer and nine-site radius gates are \(2+\sqrt{D_f-2p}<\sqrt{D_f}\) and \(\sqrt{64/9}+\sqrt{D_f-2q}<\sqrt{D_f}\). Thus the count is at most \(6+1/2+4/9=125/18\).

If no third shear \(c\) works, the forbidden balls cover the torus. Lift to Euclidean space and use the displayed replacements, intersecting the corresponding finite-index period sublattices when more than one type appears. The resulting ordinary equal-ball cover has the stated count per \(L\)-cell; coincident centres can only lower it. The exact checks \((15/2)D_g^{46}<A^{46}\) and \((125/18)D_f^{46}<A^{46}\) contradict the Coxeter–Few–Rogers threshold. Thus there is an admissible \(c\) in either branch, and the resulting lattice has squared minimum exactly 16 because \(L\) supplies vectors of that norm.

## Exact density and replay

The exact squared tail covolumes are
\[
\det G_g=\frac{35649962490781510227209447}{500000000000000000000000},
\qquad
\det G_f=\frac{14260143213204819}{200000000000000}.
\]
The second is larger (about \(71.3007161\)); both are strictly below the [pinned predecessor's worst determinant](https://github.com/alexlegeartis/SpherePackings/blob/abc44e6ad4c414b648a2d091a60f6ee08ce2f0aa/research/mw92_five_cluster95_20261128/result.json), \(3578813470795833160453541/50000000000000000000000\). The lattice covolume is \(47^7\sqrt{\det G}\), and its radius-2 occupied density is \(2^{143}\pi^{47}/(95!!\,47^7\sqrt{\det G})\). The lower decimal above is rounded outward by rational Machin intervals. The table's older reference, Shioda MW92 \(\times A_3\) of squared tail covolume \(2048\), is a derived comparison; we do not attribute a published 95D numerical lower bound to Shioda.

Run `node verify.js` (no dependencies). It replays all field equations and height pairings, the mod-two and mod-three partition criteria, both complete tail enumerations, determinant ordering, exact square-root inequalities and the decimal floor. The [earlier analytic checker](../mw92_paired_cover_9495_20261126/verify.py) independently verifies the unchanged Coxeter–Few–Rogers volume threshold \(A\) and the fixed 94D input. No optimizer or finite 95D Gram matrix is needed.
