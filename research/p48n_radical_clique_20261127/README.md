# A forbidden six-clique in one 44D radical enlargement of P48n

Sphere 143, 27 November 2026. **No new density record.** All density constants here are center densities. The live comparator in [SpherePackings main at 7cb8ce24](https://github.com/alexlegeartis/SpherePackings/blob/7cb8ce247b53b0715889178645924eb720f488cb/LOWER_BOUNDS_1_96.md) is the 44D construction of [Chen–Hu–Li–Wang–Wu](https://arxiv.org/html/2505.02394v1#S3.SS7),
\[
 \delta_0=\frac{5}{\sqrt{605}}(157/110)^{22}=509.61898264\ldots .
\]
Its occupied-volume fraction is $3.92442653352771\cdot10^{-8}$. Our two-axis canonical starting point is not the same construction.

## The new target: improve by enlarging only two radical directions

Let $P$ be the [catalogue $P_{48n}$](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ48a.html) with squared minimum 6. The two orthogonal $\sqrt 2\Lambda_{24}$ eigensublattices occur in [Nebe, Table 2](https://arxiv.org/pdf/1312.3780). The two explicit integral coefficient bases in the witness were found from the involution word $(3,4,3,4,3,4)$ in zero-based catalogue generator order, then Hermite/LLL reduced. Replay does not trust the word or LLL: it verifies the public parent Gram, two orthogonal rank-24 bases whose Grams divided by 2 are even unimodular, and their exact common mod-2 row space in the parent. Each stock is primitive: their orthogonal product has index $2^{24}$ in the parent, so no odd prime can enter a saturation index; each basis has rank 24 modulo 2, so no factor of 2 can enter either. Each stock divided by $\sqrt2$ is even unimodular with no roots because the parent has minimum 6, hence it is a Leech lattice.

Normal planes are spanned, on each side, by two stock vectors of squared length 8 and mutual dot $\pm4$. Dividing by $\sqrt2$ gives a $\sqrt2 A_2$ of Gram determinant 12; it is primitive in Leech because a proper rank-two overlattice would have determinant at most 3, below the planar Hermite bound at minimum 4. Thus the unscaled perpendicular complements $L_\pm$ are copies of $\Lambda_{22}$ of determinant 12. Use the physical period $H_+\perp H_-$ with $H_\pm=\sqrt2 L_\pm$; its covolume is $12\cdot2^{22}$.

In the common 24-bit parent parity space, each normal plane is a totally singular two-space $R_\pm$, and $H_\pm/2H_\pm$ maps onto $R_\pm^\perp$ for the Leech bilinear form modulo 2. This follows from primitive embedding, orthogonality and dimensions. The exact literal example has a nonsingular cross pairing between $R_+$ and $R_-$. In particular the two spaces are disjoint. Set $U=R_+^\perp\cap R_-^\perp$, dimension 20; the pairing condition gives $R_\pm^\perp=U\oplus R_\pm$. Phases are pairs of half-period labels. The old canonical graph has both labels equal to $u\in U$; every corresponding half sum is in the parent because the two full 48-coordinate parity rows agree. The **full radical palette** consists of
\[
    (u+r,u+s)\quad(u\in U,\ r\in R_+,\ s\in R_-),
\]
so contains $2^{24}$ distinct periods. No additive assumption is made on a selected subset of this palette.

A nonzero radical label has numerator squared length at least 12 in $L_\pm$: choose a normal root $r$ of squared length 4 in its Leech class. Its shortest perpendicular numerator $x$ obeys $(x+r)/2\in\Lambda_{24}$, hence $(x^2+4)/4\ge4$. This handles pure radical moves. For different canonical labels, neither reduced coordinate is zero; equal quadratic parity means that a physical distance below 6 could only arise from squared-4 numerators in both unscaled factors. Thus for any pair of such normal planes there are just 15 possible nonzero radical offsets to inspect. Defects are checked by matching the shallow lists: each pair of orthogonal norm-8 physical stock root lines identifies an offset via exclusive-or of their labels. Offset zero has no defects by the parent minimum. If the absent offsets contain all three nonzero elements of one two-dimensional subgroup, taking its four radical banks over **all** $2^{20}$ canonical sites gives a lattice packing of squared minimum at least 6 and
\[
 \delta=(3/2)^{22}/12=623.48563689\ldots >\delta_0
\]
(+22.34%). The search needs only the two norm-8 root lists, rather than enumerating all pairwise separations among $2^{22}$ phases. This remains a conditional constructive gate, not a found candidate. Moreover a single conflicting offset kills uniform use of that subgroup even with arbitrary deletions of canonical labels: the forbidden nonzero base difference is an involutive perfect matching, forcing at most half the canonical labels.

## Literal six-clique: even nonuniform nonlinear selection cannot fix this pair at squared minimum six

For the exact normal pair in the witness, put an edge on the full radical phase group if that difference admits a horizontal representative of squared norm 4. The checker certifies six explicit phases and all 15 edges. Each edge witness supplies one norm-8 vector from each physical stock, orthogonal to its two normals, whose half sum has squared norm $(8+8)/4=4$ and whose coefficient parities are the respective required differences. All translates of these six phases in the finite palette are thus forbidden cliques, independent of any further short-vector census.

At actual squared minimum $\mu>4$, any valid set $S$ of phases meets every translate of this clique in at most one point. Count occurrences of the phase vertices across all translates: each phase occurs six times, giving $|S|\le2^{24}/6$. In particular for $4<\mu\le6$ its center density is at most
\[
 \frac{2^{24}}{6}\frac{(6/4)^{22}}{12\cdot2^{22}}
       =\frac{(3/2)^{22}}{18}=415.65709126\ldots<\delta_0 .
\]
For $\mu\le4$, even the whole palette gives density at most $2^{24}/(12\cdot2^{22})=4/3$, also below the target. All phase differences have even squared norm: the two reduced quadratic values agree, and the physical squared length is half the sum of the even unscaled lengths. Thus the only higher-minimum possibility is $\mu=8$, since the period contains norm-8 vectors. The theorem controls arbitrary nonlinear motifs on this **one literal pair of normals**, fixed period/scale, and squared minimum at most 6; it does not cover $\mu>6$, perturbing the geometry, non-binary phases, or other normal pairs.

Exploration found this clique after exact LDL enumeration of both eigensublattice shortest shells (196,560 vectors per stock, with one zero shorter in each) and a seeded $20\times20$ plane sample. Only 155 of these pairs had nonsingular cross form. Every nonzero offset had a collision in each of those 155, but this is not an exhaustive normal-orbit claim. The literal proof uses **only 15** exact short representatives, not the exploratory census.

## Alternative construction routes examined

* In a genuinely different ternary polarization, $L=\Lambda_{22}$ has mod-3 radical of dimension one. A normal representative is $r_1-r_2$ of squared length 12 where the two squared-4 roots have dot $-2$. Any perpendicular numerator in this nonzero radical class has squared length at least 24, because its sum with the normal divided by 3 is Leech. Thus all nine radical choices in two factors are safe, and a rank-20 anti-isometry on the 21-dimensional nondegenerate quotients for two $\sqrt3 L$ factors would reach the same $\delta=623.4856\ldots$ if its three shallow combinations $(4,8),(6,6),(8,4)$ are absent. Off-diagonal norms cancel modulo six, and one-sided nonzero radical classes remain safe. A structured source for the mod-3 match is [Coulangeon–Nebe's Hermitian $P_{48n}$ tensor](https://arxiv.org/pdf/1201.1832), with $\eta=(1+\sqrt{-11})/2$. In its binary Hermitian matrix $\left(\begin{smallmatrix}2&\eta\\\bar\eta&2\end{smallmatrix}\right)$, $(1,-1)$ and $(2-\eta,2-\bar\eta)$ are orthogonal norm-3 vectors (using conjugate-linearity on the first variable), and their determinant is the scalar 3, of additive index 9. In that frame the full lattice uses the relation $x=(1+\eta)y\bmod3$. This section-level lifting still needs explicit shallow separation. The exact Golay-coordinate theta scout counted 49,896, 2,821,632, and 50,332,590 points at squared norms 4, 6, 8 in $\Lambda_{22}$; a random anti-isometry is not supported by these counts.
* Independent reserve: release one primitive dual direction from [Dorofeev–Sun–Wang's 43D $P_{48m}$ normal](https://arxiv.org/html/2607.20359v4#S4). In the finite range $48\|w\|^2\le48$ for the released covector and normal squared diameter $\beta<2/3$, 130 sign-normalized primitive releases yield no improved rank-four anticode certificate. The best seen is the already-known five-site $7/24$ diameter over the determinant-768 section, center density $450.99813389\ldots$. This cannot rule out stronger actual horizontal radii even for those periods; outside the stated finite range it is not a bound. Neither exploration produced an improved packing.

Replay the literal theorem with `node verify.js`. The mathematical minimum of the public parent is a literature input. All comparisons and 15 edges are exact integers; no numerical norm estimate enters the proof.

