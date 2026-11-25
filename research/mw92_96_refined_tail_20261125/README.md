# An exact metric/color refinement of the 96D twelve-color MW92 packing

This is a periodic existence construction with occupied density
\[
\delta_{\rm new}=\frac{\pi^{48}2^{96}}
 {48!\,47^7\sqrt{300153807922100058481/1440000000000000000}}
 \ \ge\ 6.367793872884866411\cdot10^{-22}.
\]
It improves the preceding verified **unrestricted** lower bound by approximately **0.00030269067%**. The change is deliberately reported as a small exact metric refinement, not a new horizontal construction, a lattice packing, an upper bound, or a claim about literature priority.

## Inherited horizontal theorem and comparator

The [preceding exact construction](../mw92_96_periodic12_20261124/) proves the rank-92 MW carrier \(L\) has covolume \(47^7\) and minimum squared 16. Its certified commuting order-three isometry \(J\) supplies \(H=\{0,h,Jh,J^2h\}\) as a subgroup modulo \(L\), with every nonzero element at squared torus depth at least \(21/2\). The independently certified four-ball covering argument yields \(b\) for which \(b+H\) is at squared depth \(>106743/10000\). The argument in the preceding construction has already verified that the twelve offsets \(H+\{0,b,-Jb\}\) have that same bound across different ternary groups, by applying \(J\) and sign. These statements, including the entire half-phase CRT enumeration and the four-ball inequality, are reused unchanged. The [pinned comparator result](https://github.com/alexlegeartis/SpherePackings/blob/50c5bc5f66c88a8625bfd74b4a045c3db1bf2752/research/mw92_96_periodic12_20261124/result.json) has squared tail covolume \(53361/256\) and occupied density floor \(6.367774598225213\cdot10^{-22}\). Its exact density is \(\pi^{48}2^{96}/(48!47^7\sqrt{53361/256})\).

## New tail and coloring

In ordinary integer coordinates \(v=(i,j,k,l)\) use \(D_4=\{v:\ i+j+k+l\equiv0\bmod 2\}\). Its basis is
\[
(1,-1,0,0),\quad(0,1,-1,0),\quad(0,0,1,-1),\quad(0,0,1,1).
\]
Write \(v\) in that basis using coefficients \((a,b_0,c,d)=(i,i+j,(i+j+k-l)/2,(i+j+k+l)/2)\). The new color is
\[
(c+d\bmod 3,\ a\bmod 2,\ b_0+c\bmod 2).
\]
The map is onto by independent mod-3 and mod-2 rank, or by the explicit residue replay. Identify the binary group with \(H\); to its element add the ternary offset \(0,b,-Jb\) as indicated. Use the following Gram on **ordinary** \((i,j,k,l)\) coordinates, not basis coefficients:
\[
Q=\frac1{120000}
\begin{pmatrix}
325000&-4771&-5229&-458\\
-4771&325458&-5000&229\\
-5229&-5000&324542&229\\
-458&229&229&315000
\end{pmatrix}.
\]

At equal ternary but different binary color the vertical target is \(16-21/2=11/2\); at distinct ternary color it is \(16-106743/10000=53257/10000\). At the same color the target is 16. Gershgorin applied to the integer matrix gives \(\lambda_{\min}(Q)\ge314084/120000>16/7\). Therefore every nonzero tail difference of ordinary norm squared at least 7 automatically meets all three thresholds. Each tail difference of ordinary norm squared at most 6 has every coordinate in \([-2,2]\), giving a complete small integer check:

| Difference color | Number of signed tail vectors with ordinary squared norm \(\le6\) | Minimum numerator, denominator 120000 |
| --- | ---: | ---: |
| distinct ternary | 90 | 639084 |
| same ternary, distinct binary | 42 | 660000 |
| same color | 12 | 1920000 |

Every distinct vertical fiber separation thus adds up with its corresponding inherited horizontal separation to at least 16. Within a fiber the horizontal minimum 16 is attained in the carrier. The result uses twelve cosets over \(L\) times the kernel of the tail color; the kernel index 12 exactly cancels the motif cardinality. The ordinary coordinate \(D_4\) lattice has covolume 2 before the Gram transform; its new squared covolume is
\[
4\det Q=\frac{300153807922100058481}{1440000000000000000}<\frac{53361}{256}.
\]
Taking radius 2 gives the stated density. The exact comparison does not use the displayed decimal: \(256\cdot300153807922100058481<53361\cdot1440000000000000000\).

## Replay

Run `node verify.js`. It verifies complete tail coverage, the onto color map, every constraint and minimum numerator, Gershgorin for long differences, strict exact determinant comparison, and the occupied-density lower floor using a Machin-series lower bound on \(\pi\). The prior certified horizontal proof is linked above and is not recomputed by this short tail-only script. The lattice-only lower table remains unchanged.
