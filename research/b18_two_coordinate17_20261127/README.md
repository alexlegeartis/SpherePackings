# Two-coordinate shadows of the Bierbrauer–Edel packing in 17D

Sphere146, 27 November 2026. **No new density record.** The benchmark is the laminated \(\Lambda_{17}\) with center density \(1/16\), occupied fraction \(0.00881131918232118\ldots\); see [Leech–Sloane, §4.5](https://neilsloane.com/doc/Me14.pdf) and the [live project comparison](https://github.com/alexlegeartis/SpherePackings/blob/main/LOWER_BOUNDS_1_96.md). The parent construction \(B_{18}\), not ours, is due to [Bierbrauer–Edel, §2](https://www.yvesedel.de/Papers/Sphere.pdf) and has center density \((3/4)^9\). These results concern specific 17D projections, *not* arbitrary 17D lattices or packings.

Two distinct gains were plausible: use both heights of every projected parent phase (doubling the naive shadow), or free the real axial positions and high-digit parities of every exterior QR word. Both are excluded at **all** minima by different arguments. The first gives a literal 1,024-site packing at 98.0902% of the benchmark; the second proves, even with arbitrary real axial phases, an all-minimum family ceiling below 99.532%. This leaves changes to the exterior metric or alphabet, not another scalar or axial search.

## 1. Parent, projection, and exact geometry

Use both extended QR codes \(C_0,C_1\) of length \(18\): form the cyclic binary spans of the quadratic-residue and nonresidue masks on \(\mathbb Z/17\), append a last zero bit, and adjoin the all-ones 18-bit word. Each is a 512-word even linear code of minimum Hamming distance six; they are orthogonal to each other. With \(P=4D_{18}\), the two parent banks are \(2a+P\) for \(a\in C_0\), and \(\mathbf 1+2b+4e_o+P\) for \(b\in C_1\), where \(e_o\) is any coordinate vector. Indeed same-bank squared distance is \(4d_H\) for distinct words (one \(\pm2\) roundoff repairs the even-sum condition of \(D_{18}\)), and cross-bank pairs must pay one \(3^2\) instead of \(1^2\) because the required high-digit parity is \(1+|b|+a\cdot b=1\bmod 2\). Thus their squared distance is 26. Period minimum is 32; center density is \(1024\cdot 6^9/(2\cdot 4^{18})=(3/4)^9\).

Delete the normal direction \(e_i\pm e_j\), keep the sixteen exterior coordinates, and write the remaining axial coordinate as \(x_i\mp x_j\). Use orthogonal metric \(t/2\) times squared axial coordinate, \(t>0\), and ordinary raw Euclidean metric on the exterior. The period is the horizontal intersection \(P_0\), with separate axial period 8 and exterior period \(4D_{16}\). This holds for either sign because the two integral active coefficients in \(P\cap(e_i\pm e_j)^\perp\) have even sum. Its covolume is \(2\cdot4^{17}\sqrt{2t}\). The projection of \(P\) itself has exactly two cosets modulo \(P_0\): a generator of the other is a vector with axial increment \(4\bmod 8\) and exterior odd high digit. Indeed normal pairings of \(P\) are the multiples of 4, and the primitive vertical vector in \(P\) has normal coordinate 8. Thus the full stock has 2,048 phase sites, two per parent word; distinct words in the same bank have distinct exterior binary residues because the parent minimum is six.

Each phase lies in one of four raw axial labels \(0,2,4,6\bmod 8\). Equal, adjacent and opposite labels contribute respectively \(0,2t,8t\) to squared torus distance. In either bank, the parity of the label divided by 2 agrees with the sum of the two active bits, up to a bank-dependent constant. For same-bank distinct words whose exterior Hamming distance is \(d\), the exterior squared distance is exactly \(4d\) regardless of high-digit parities: every nonzero \(d\) includes a \(\pm2\) ambiguity. For the two lifts of the same word it is \(16\). Cross-bank exterior squared distance is \(16\) or \(24\), according to whether the nearest signed unit residues have the prescribed parity. All these are torus distances on a product, so the two parts minimize independently. The period enforces \(\mu\le \min(32,32t)\) for the actual squared minimum \(\mu\). If \(M\) sites are chosen, the exact density ratio to \(1/16\) is
\[
R=\frac{M}{64\sqrt t}\left(\frac{\mu}{32}\right)^{17/2}.       \tag{1}
\]

**Rectangle fact.** In either code and every choice of the two active coordinates there are weight-six words \(g_1,g_2\) with active projections \(10,01\), such that \(g_1+g_2\) also has weight six. This finite fact checks all \(2\binom{18}{2}=306\) choices directly in [verify.js](verify.js). The cosets of their two-dimensional span partition the bank into 128 four-word rectangles. The four side differences have exterior Hamming weight five, and the two diagonals have exterior Hamming weight four. In particular \(g_1+g_2\) also supplies a matching between the two active patterns of each equal label parity.

## 2. The all-metric theorem for the full stock

**Theorem.** For every \(t>0\), every subset of the 2,048 projected phases modulo \(P_0\), at its actual minimum, has \(R\le1\). This includes arbitrary nonlinear and partial-bank subsets; no assumption of one lift per word or full binary fibers is made.

If \(\mu\le16\), using \(M\le2048\) and \(t\ge\mu/32\) in (1) gives \(R\le1/8\). If \(16<\mu\le24\), each fixed axial label and bank admits at most 128 sites: pair its 256 possible words by \(g_1+g_2\); matched torus squared distance is 16. Thus \(M\le1024\). When \(\mu\le20+2t\), insert \(\min(32t,20+2t,24)\) for \(\mu\) in (1). The ratio bound increases until \(t=2\), then decreases: the switch from \(32t\) to \(20+2t\) is at \(2/3\), and \((20+2t)^{17/2}/\sqrt t\) increases for \(t>5/8\). At \(t=2,\mu=24\) the ratio is
\[
8\sqrt2\,(3/4)^{17/2}=\sqrt{128(3/4)^{17}}=0.9809022340\ldots .
\]
If instead \(\mu>20+2t\), all four rectangle sides are forbidden even under arbitrary choices of lifts, since the ends have adjacent labels. Consequently either only one word appears, or only the two diagonal words. Each diagonal pair permits at most two sites in total, by matching its identical labels (distance 16). This gives \(M\le512\) and \(R\le8(3/4)^8<1\).

For \(\mu>24\), each fixed label may contain only one bank (mixed exterior distance at most 24). Within it the length-16 exterior binary code has minimum at least 8: distances are even because active parity is fixed, and exterior squared distance is \(4d\). Therefore \(M\le4A(16,8)\le128\). The elementary bound \(A(16,8)\le32\) follows by mapping to unit sign vectors with nonpositive pairwise inner products; the sum of squares of off-diagonal Gram entries is at most their absolute sum, at most the trace, and trace–rank then gives \(M\le 2\cdot16\) in one label.

* If \(t\ge4\), (1) is at most 1 using \(M\le128,\mu\le32\).
* Suppose \(t<4\). Then \(t>3/4\) by \(\mu>24\) and the period bound. If \(\mu\le24+2t\), (1) is at most its value for \(M=128,\mu=24+2t,t=4\), namely 1: the enveloping expression is increasing for \(t>3/4\).
* Otherwise \(\mu>24+2t\) and a mixed-bank pair can use **only opposite** labels. If both banks appear, each occupies a single label, so \(M\le64\) and \(R\le1\) by \(\mu\le\min(32,32t)\).
* With a single bank, a union on each adjacent pair of labels is a length-16 code of minimum at least 7: an adjacent-label pair has axial contribution \(2t\), and \(\mu>24+2t>20+2t\). There are two disjoint adjacent pairs, so \(M\le2A(16,7)=72\), by [Östergård's 2011 result](https://doi.org/10.1109/TIT.2011.2144955). For \(t\ge81/64\), this suffices directly in (1).
* For \(3/4<t<81/64\) and \(\mu\le20+8t\), the bound with \(M=72,\mu=20+8t\) increases with \(t\); its value at \(81/64\) is smaller than 1 because \(20+8(81/64)<32\). Otherwise \(\mu>20+8t\). Append two Gray bits for the four axial labels, so equal, adjacent, opposite labels have extra Hamming distances \(0,1,2\). Exterior distances are respectively \(\ge8,\ge7,\ge6\); in the opposite case the axial squared contribution is \(8t\), and even repeating a word in both lifts would cost just \(16+8t\), already forbidden. This produces a length-18 binary code of minimum 8. [Östergård's \(A(18,8)=64\)](https://acris.aalto.fi/ws/portalfiles/portal/27672730/ELEC_ostergord_golay.pdf) gives \(M\le64\), sufficient as before.

For a literal near-benchmark construction choose the plus normal at the first and last coordinates, \(t=2\), and number the active patterns \(h=w_i+2w_j\). Starting at axial label \(2(w_i-w_j)\) and a fixed odd exterior high digit for the odd bank, lift the even bank by \(4\) plus an odd exterior high digit exactly when \(h=0\), and the odd bank exactly when \(h=3\). Keep all 1,024 resulting words. The checker gives minimum squared distance 24, with 112,640 unordered motif pairs attaining their torus minimum at 24. This is a new literal construction in the stated family, not a priority, isometry, or density-record claim.

## 3. Free axial phase, period, and high digit, through \(\mu=24\)

Discard the four axial labels and the lift rule entirely. Retain the fixed exterior period \(4D_{16}\) and the projected binary words of the two banks. Allow any axial period \(T\), any number of sites per word, any real axial positions, and either exterior high-digit parity independently for each site. Let \(\rho=M/T\). The density ratio is now
\[
R=\frac{\rho}{8\sqrt2}\left(\frac{\mu}{32}\right)^{17/2}.
\]
For any fixed bank, word and high parity, axial positions must be \(\sqrt{\mu}\)-separated modulo \(T\) (and the period itself must be at least \(\sqrt{\mu}\)). Hence \(\rho\le2048/\sqrt{\mu}\). In each of the 256 rectangles across both banks, any two selected sites have same-bank exterior squared distance at most 20: for different words use the rectangle differences; for the same word it is 0 or 16. Thus for \(\mu>20\), \(\rho\le256/\sqrt{\mu-20}\). No cross-bank constraints are needed. For \(\mu\le41/2\) the first bound gives \(R^2\le 1024(41/64)^{16}<1\). On \([41/2,24]\) the second bound gives \(R^2\le512(\mu/32)^{17}/(\mu-20)\); its only interior critical point (\(\mu=85/4\)) is a minimum. The endpoint bounds \(1024(41/64)^{17}<1\) and \(128(3/4)^{17}<1\) certify the entire remaining interval. In particular no coordinated compression of the displayed 1,024-site pattern, even with unrestricted axial positions, high digits, deletions and multiples, improves it at **the same raw squared minimum 24**.

## 4. Finite window code bounds close the rest of the free-axis family

**Stronger theorem.** The free real axial stock in §3, with arbitrary period, multiplicities, real positions and high-digit parities, in fact satisfies at every minimum
\[
R^2\le126(379/504)^{17}<1,\qquad R\le0.995311443826\ldots .
\]
The same upper bound holds for upper asymptotic density without an axial period, with the same complete exterior translates. The exterior period forces \(\mu\le32\).

**Why the punctured QR stock is tighter than arbitrary binary codes.** For either bank \(C\), let \(C'\) be the length-16 puncture retaining the exterior coordinates. For every choice of the two deleted coordinates, the largest subsets of \(C'\) of Hamming minimum 7 and 8 have sizes **32 and 16**, respectively. This is a finite exact clique result. First, cyclic shifts and the coordinate map \(x\mapsto-1/x\) on \(\mathbb F_{17}\cup\{\infty\}\) preserve both length-18 codes (readily verified); they move any ordered pair to \((0,\infty)\): first send one position to \(\infty\), then translate the other to zero. It suffices to puncture that pair. The puncture remains a binary linear code of size 512, so any nonempty subcode can be translated by a member to contain zero. Filter the other codewords to those whose weight is at least the requested minimum, and join words whose difference has sufficient weight. A maximum clique among these neighbors verifies each bound. The explicit maximum witnesses and the branch-and-colour verifier are in [verify.js](verify.js); the colouring is only an exact branch upper bound, not a numerical clique heuristic. For minima 7 and 8 the upper search excludes cliques of size 32 and 16, respectively, among the nonzero neighbors, and the supplied neighborhoods of size 31 and 15 certify equality.

**Apply to short axial intervals, not just fixed layers.** An interval of length \(\sqrt{\mu-24}\) when \(\mu>24\) can be taken half-open. Any two of its axial positions then differ by strictly less than this length, requiring exterior squared separation greater than 24. The interval uses only one bank (cross-bank exterior gap never exceeds 24), no repeated exterior word, and is a punctured-QR code of Hamming minimum at least 7. It therefore contains at most 32 sites. Similarly when \(\mu>28\), an interval of length \(\sqrt{\mu-28}\) contains at most 16 sites (exterior Hamming minimum at least 8). Divide an arbitrarily long interval into such windows to get
\[
\rho\le \frac{32}{\sqrt{\mu-24}}\quad(\mu>24),\qquad
\rho\le \frac{16}{\sqrt{\mu-28}}\quad(\mu>28),
\]
with uniform count errors independent of the initial location. The elementary bounds of §3 also apply without an axial period by counting per-label spacings in arbitrary intervals. Each corresponding full horizontal translate has the same density and a uniform geometric-boundary count error; summing over slabs, or approximating an expanding ball by horizontal boxes and slabs, carries the uniform bound on \(\rho\) to the stated spatial upper density.

For the exact all-\(\mu\) computation use the four bounds for \(R^2\) on the following four closed intervals, assigning shared endpoints to either:
\[
\begin{array}{c|c}
 (0,41/2]&1024(\mu/32)^{16}\\
 [41/2,1516/63]&512(\mu/32)^{17}/(\mu-20)\\
 [1516/63,88/3]&8(\mu/32)^{17}/(\mu-24)\\
 [88/3,32]&2(\mu/32)^{17}/(\mu-28).
\end{array}
\]
The middle boundaries are the respective crossings. All expressions of the form \(c(\mu/32)^{17}/(\mu-a)\) have only the critical point \(17a/16\), a minimum, so taking endpoints suffices. The largest endpoint value is at \(1516/63\), giving \(126(379/504)^{17}\). The exact comparisons are included in the checker. In particular at \(\mu=32\) the separate bound is \(R\le1/\sqrt2\). An improvement in this direction now needs to change **exterior** stock or metric, not the axial spacing schedule.

## 5. Two independently weighted exterior orbits in the literal pattern

For the displayed 1,024 sites, keep all point coordinates and split the exterior coordinates into the eight quadratic residues and eight nonresidues. Give the two orbits independent squared metric factors \(a,b>0\), and give the physical axial coordinate its own squared factor \(c>0\) (so the earlier \(t=2\) becomes \(t=2c\)). Even this distinct deformation cannot improve the literal pattern. Normalize its minimum to at least 24 by common scaling. The shortest squared motif distances are linear forms in \(a,b,c\), except a cross-bank odd parity repair uses \(8\min(a,b)\), which is equivalent to imposing both linear inequalities. After dropping dominated forms the complete period/motif conditions are the twelve inequalities whose coefficient rows are
\[
(32,0,0),(0,32,0),(0,0,64),(8,16,0),(16,8,0),
(4,16,4),(16,4,4),(12,8,4),(8,12,4),
(8,8,16),(16,0,16),(0,16,16),
\]
each with lower bound 24. The complete vertices are \((1,1,1)\), \((3/2,3/4,3/2)\) and its \(a,b\) swap, \((9/8,9/8,3/8)\), and \((3/4,9/4,3/4)\) and its \(a,b\) swap. The recession cone is the nonnegative orthant. The relative covolume is \((ab)^4\sqrt c\); its logarithm is concave and coordinatewise increasing, so checking the listed vertices suffices. The minimum is 1. The exact checker recomputes the whole motif pair support and the vertices with integer arithmetic. This does not treat arbitrary exterior metrics or changing sites there.

Replay of the exact finite ingredients: `node verify.js`. The two published binary-code upper bounds used in Section 2 are mathematical inputs, not re-proved by the script.
