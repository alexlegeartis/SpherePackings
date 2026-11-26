# 17D: the Hamming–\(D_4\) deep-hole stagger gate

November 26, 2026. This is an exact bound for a particular stacking mechanism, with an exact equality construction, **not a sphere-packing record**. The pinned [project table](https://github.com/alexlegeartis/SpherePackings/blob/8b652d31db18077dc97415aa17d87f52ffba01a0/LOWER_BOUNDS_1_96.md) and [Cohn's table](https://cohn.mit.edu/sphere-packing/#data) give center density \(1/16\) for 17D \(\Lambda_{17}\) (ordinary occupied fraction approximately \(0.00881131918232118\)). Use minimum squared distance \(4\) throughout, so center density is simply the number of centers per unit volume. For historical code-based layering context see [Conway and Sloane, What are all the best sphere packings in low dimensions?](https://neilsloane.com/doc/BEST/Me194.pdf); no priority or isometry claim is made for the realization below.

## The 16D sheet

Partition 12 coordinates into three blocks of 4, with positions in each indexed by \(\mathbb F_2^2\). The binary check map keeps the three block parities and the \(\mathbb F_2^2\) sum of the selected positions across all blocks. Let \(C\) be its kernel and \(H=C+2\mathbb Z^{12}\). The 12 columns are distinct, nonzero, and a sum of three has at least one block-parity bit nonzero. They span the five check bits. Thus \(C\) has size \(128\), \(H\) has minimum squared distance \(4\) (attained by \(2e_j\)) and covolume \(32\).

A **chart phase** for a translate of \(H\) assigns each block a symbol \((\alpha,p)\in\mathbb F_2^2\), meaning a shift by \(\alpha/2\) in all four coordinates and integer syndrome block parity \(p\); a remaining tag \(t\in\mathbb F_2^2\) specifies the total two-bit position sum. A doubled half-shift in one block has zero syndrome, so the phase parameters add over \(\mathbb F_2\). If two chart phases have distinct symbols in \(k>0\) blocks, their squared separation is \(k\): a block with unmatched half-offset costs \(1\), with either sign on each half-coordinate, so can absorb arbitrary tag changes; a block with equal half-offsets but different prescribed parity also costs \(1\), and the position of a single needed odd coordinate can absorb any tag change. If all symbols coincide, the squared separation is \(0\) or \(2\) according as the tags coincide or differ. In the latter case two distinct positions in one block realize any nonzero tag.

A full sheet places a copy of \(H\) with a chart phase at every normal site in \(D_4/\sqrt 2\), with \(D_4=\{v\in\mathbb Z^4:\sum v_j=0\bmod 2\}\). Allow any phase assignment as long as the sheet has minimum squared distance at least 4; neither additivity nor a particular period is imposed. The four-dimensional normal grid has covolume \(1/2\), so the sheet's center density is \(1/(32/2)=1/16\). One such legal sheet is exhibited below.

## The arbitrary-horizontal-shift bound at a deep-hole stagger

A standard deep-hole stagger in \(D_4/\sqrt 2\) has squared depth \(1/2\), with eight nearest sites forming a cross-polytope. In integer coordinates before scaling these staggers are an odd-sum integer vector or an all-half-integer vector. In the integer case the eight nearest points are the four pairs \(v\pm e_i\); in the half-integer case they are the qualifying eight binary corners of a unit cube, paired by complementation. In both cases opposite vertices have normal squared distance \(2\), and other distinct pairs have normal squared distance \(1\). The graph joining the latter pairs is \(K_{2,2,2,2}\). To reach packing separation \(4\) on these edges the horizontal separation must be at least \(3\); the chart metric forces all three block symbols to differ. Consequently each block uses its four symbols once each on the four opposite pairs, with a common symbol within each pair. The two inner tags in an opposite pair must differ (the normal cost is only \(2\)).

**Claim.** An arbitrary point \(g\in\mathbb R^{12}\) is horizontally at squared distance at most \(5/2\) from one of the eight translates above. In fact the conclusion uses only the four distinct block symbols and the distinct opposite tags just derived.

For each opposite pair \(i\), the two allowed tags form a coset of a line in \(\mathbb F_2^2\), so their union specifies only one linear tag bit. This bit partitions each block's four positions into two coordinate pairs, and the partition may be different for each \(i\). In one block write \(u_j=2\,\operatorname{dist}(g_j,\mathbb Z)\in[0,1]\). Under the two half-offsets the normalized residuals are \(u_j(0)=u_j\) and \(u_j(1)=1-u_j\). Fix \(\alpha\); rounding in that shift gives base squared cost \(S_\alpha=\sum_j u_j(\alpha)^2/4\). Flipping the rounded parity in the two coordinate pairs costs respectively \(a=1-\max_{j\text{ in first pair}}u_j(\alpha)\) and \(b=1-\max_{j\text{ in second pair}}u_j(\alpha)\). These are exact: the least loss is a single parity change at the coordinate with largest residual. When the prescribed block parity agrees with rounding, costs for the two retained tag bits are \(S_\alpha\) and \(S_\alpha+a+b\) in some order; if it disagrees, the costs are \(S_\alpha+a\) and \(S_\alpha+b\). Denote the lower cost by \(l\) and the gap by \(g'\). The global retained bit is the XOR of the three block contributions, so the nearest translate in pair \(i\) costs at most
\[
F_i=\sum_{b=1}^3 l_{ib}+\min_b g'_{ib}.
\]

For a block set \(m=\min u_j\), \(M=\max u_j\), \(R=M-m\), \(E=\sum_j u_j(1-u_j)\), and \(D=R+E\). Each of the four \((\alpha,p)\) labels appears once, giving
\[
\sum_i l_{ib}=2S_0+2S_1+2-M-(1-m)=3-D.
\]
This does not depend on any of the pair partitions: the lower loss for the disagreeing parity is always \(1-\max_j u_j(\alpha)\). The disagreeing-parity gaps are at most \(R\). The agreeing-parity gaps for the two half-offsets are at most \(2-m-M\) and \(m+M\), respectively. The smaller of these bounds is at most \(D\): for \(m+M\le1\), every \(u_j\in[m,1-m]\), giving \(E\ge 4m(1-m)\ge 2m\), hence \(R+E\ge m+M\). Reflect \(u\mapsto1-u\) for the other case. Thus at most one label in each block can have gap exceeding \(D\).

If some \(i\) has \(d=\min_b g'_{ib}>1\), it must agree with rounding in every block. In an agreeing block with pair maxima \(U,V\), \(l=S_\alpha\le(U+V)/2=1-g'/2\), by \(u^2\le u\) and two coordinates per pair. This gives \(F_i\le3-3d/2+d=3-d/2<5/2\). Otherwise every pair's minimum gap is at most \(1\). Select the block of least \(D\); there at most one exceptional pair may have gap above \(D\). Hence
\[
\sum_i F_i\ \le\ 9-\sum_b D_b + 1+3\min_b D_b\ \le\ 10.
\]
At least one \(F_i\le5/2\), proving the claim, including arbitrary shifts \(g\) and ties in rounding.

For any upper sheet in a layer stack, choose one center above a normal deep hole of the adjacent lower sheet. The 16D projected squared distance from some lower center is at most \(1/2+5/2=3\); radius-one packing forces last-coordinate gap at least \(1\). This is pointwise and allows arbitrary real horizontal shifts between sheets, unrelated nonlinear phase assignments, arbitrary actual last-coordinate positions provided consecutive sheets have one of the deep-hole staggers, and no common period. Counting sheets (at most one per unit length asymptotically) gives 17D center density at most \(1/16\). For nonperiodic assignments, use the uniform counting error of translates of fixed \(H\) and \(D_4/\sqrt2\), or count in expanding coordinate boxes first. The bound does not apply to arbitrary stagger vectors or to sheets that alter the common coordinate chart, delete grid sites, or replace copies of \(H\) by smaller motifs.

## Exact equality construction

This is additive: the normal grid is \(y\in D_5\), mapped to \((y_1/\sqrt2,\ldots,y_4/\sqrt2,y_5)\). For \(z=(y_1,y_2,y_3,y_4-y_5)\in D_4\), put at \(y\) the translate of \(H\) specified by the following phase:
- the three block symbols are equal, and the common two-bit symbol is the XOR of the indices \(j=0,1,2,3\) of odd \(z_{j+1}\) (low bit is half-shift flag);
- the tag is \((\sum z_j/2\bmod2)\mathbin{\mathrm{XOR}}2(z_4\bmod2)\);
- add \(y_5\) times the half-grid phase \(\eta=26497=32\cdot828+1\). In this representation a phase \(32a+s\) sets half-coordinate mask \(a\) and integer syndrome \(s\). The sum of phases adds carry syndrome of the intersection of their masks.

These prescriptions define a homomorphism \(D_5\to(\tfrac12\mathbb Z^{12})/H\). The translate within each fixed last-coordinate level is a chart-decorated sheet, up to the common real shift \(y_5\eta\), and consecutive normal \(D_4\) grids differ by a deep-hole stagger. The horizontal minimum squared distance for a phase \(32a+s\) is exactly \(|a|/4+\min_{\operatorname{syn}(x)=s}|x\setminus a|\): squared cost on fractionally shifted coordinates is \(1/4\) regardless of binary choice; other odd coordinates cost \(1\). The 320 nonzero \(y\in D_5\) with normal squared norm below \(4\) all satisfy
\[
2\sum_{j=1}^4 y_j^2+4y_5^2+4\,\operatorname{dist}(32a+s,H)^2\ \ge\ 16 ,
\]
where the phase is the one prescribed at \(y\). The adjacent-layer shell of norm squared \(3/2\) attains equality, showing sharpness of the \(5/2\) horizontal gate. The finite check in [verify.py](verify.py) is independent of the search or any floating-point tolerances. At \(y=0\) the horizontal minimum is already \(4\), and all other normal vectors have normal norm squared at least \(4\). The normal covolume is \(2(1/\sqrt2)^4=1/2\); the 17D covolume is \(32/2=16\).

Replay: `python3 verify.py`.

## Why the scope matters

This rejects the direct idea of squeezing equally dense 16D Hamming–\(D_4\) sheets by a horizontal shift at a nearest-layer deep-hole stagger, independent of the alignment of their phase decorations. It does not decide anisotropic off-axis staggering or irregular normal grids. A benchmark for an off-axis transfer is finding an intersheet squared separation strictly above \(3\) before inserting the seventeenth coordinate. A large hole of one sheet is necessary but not sufficient: every translate of the next sheet, not merely one, must fit.
