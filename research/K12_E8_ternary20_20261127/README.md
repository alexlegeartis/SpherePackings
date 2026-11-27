# A full ternary Coxeter–Todd palette and two 20D stacking bounds

Sphere148, 27 November 2026. **No new lower bound.** The [Chen–Hu–Li–Wang–Wu 20D antipode](https://arxiv.org/html/2505.02394v1#S3.SSx3) has exact center density \(3^{20}/(2^{10}5^{21/2})=0.1559337388704\ldots\), and the [project table](https://github.com/alexlegeartis/SpherePackings/blob/main/LOWER_BOUNDS_1_96.md) gives its occupied density \(4.02416506306819\cdot10^{-3}\). The proposed contrasting construction keeps full horizontal Coxeter–Todd fibers, assigns real translates in the **whole** ternary palette \((K_{12}/3)/K_{12}\), and takes vertical sites in \(\sqrt{s}E_8\). At \(s=2/3\) and squared minimum 4, using all sites would have given center density \(3/16\), about 20.24% above the comparator. The obstruction below prevents this, even with nonperiodic deletions and labels.

**Theorems.** Keep \(K_{12}\) at squared minimum 4, permit any subset of vertical sites, and allow multiple horizontal phases at a site as long as the packing separates.
* If \(s\ge5/8\), every such packing, with or without a vertical period, has center density at most
\[
\frac79\frac{(47/48)^{10}}{27(5/8)^4}=0.15294592227\ldots<\delta_{\mathrm{CHLWW},20}.
\]
* If all assignments repeat under \(2\sqrt{s}E_8\), without a horizontal label shift, then for **every \(s>0\)** and squared minimum \(\mu>8/3\) the center density is at most
\[
\frac{102}{432}\left(\frac89\right)^6=0.116458\ldots<\delta_{\mathrm{CHLWW},20}.
\]
The compact-family exclusion also holds at all radii for \(s\ge5/12\), with a uniform cap \(7(2/3)^{10}/[27(5/12)^4]=0.14916490372\ldots\). These statements do not exclude low-radius compact cases with smaller normal scale, arbitrary off-palette horizontal shifts, other metrics, or general 20D packings.
* At the desired **inherited squared minimum 4** no periodicity is needed all the way down to \(s\ge7/20\). If \(s<2/3\), horizontal clearance cannot reach the required root-edge distance, so the E8 projection bound gives \(D(4,7/20)/16=0.15425671402\ldots\). If \(2/3\le s<1\), the nine-point window gives at most \(7/48=0.1458333\ldots\). If \(s\ge1\), even full occupancy gives \(1/27\). All are below the benchmark.

## 1. Exact alphabet: half a million labels, but a gap before the deepest shell

Write \(\omega=e^{2\pi i/3}\), \(R=\mathbb{Z}[\omega]\) and \(|a+b\omega|^2=a^2-ab+b^2\). Use the normalized construction \(K=\{x\in R^6:x\bmod2\in H_6\}\) of \(K_{12}\); see [Conway–Sloane](https://neilsloane.com/doc/Me107.pdf), Sections 2–3. In \(R/2R=\mathbb{F}_4\), write its elements as \(0,1,\omega,\omega^2\). Our hexacode is
\[
(u,v,w,u+v+w,u+\omega v+\omega^2w,u+\omega^2v+\omega w).
\]
Its weight distribution is \(1+45z^4+18z^6\); it is Hermitian self-dual. These properties certify the minimum 4 (the vectors \(2e_i\) give equality), covolume \((\sqrt3/2)^6 4^3=27\), and \(K^*=\theta^{-1}K\) with \(\theta=1+2\omega\) and \(\theta^2=-3\). For the dual identity, \(\theta^{-1}K\) pairs integrally with \(K\) because Hermitian pairings lie in \(2R\), and the covolumes then agree. Thus dual phases are among the ternary phases.

The \(3^{12}\) ternary labels are exactly the pairs of coefficients \(a_i,b_i\bmod3\): every residue vector lifts to \(K\) by the Chinese remainder theorem, and \(K\cap3R^6=3K\). For each coordinate residue modulo 3 and each hex symbol modulo 2, minimize \(a^2-ab+b^2\) in the corresponding pair of residues modulo 6. This is a finite exact minimization: a centered representative has norm at most 27, and the norm is at least \(3a^2/4\) and \(3b^2/4\), so \(|a|,|b|\le6\) suffices. The sum for each of the 64 hexacode words gives all horizontal torus distances; divide the minimum score by 9.

The full depth histogram in raw units is:
\[
\begin{array}{c|rrrrrrrrrr}
9d^2&0&4&6&8&10&12&14&16&18&24\\ \hline
\#&1&756&4032&20412&60480&119196&163296&122472&40544&252 .
\end{array}
\]
In particular no ternary phase has squared distance in \((2,8/3)\). **Every** label of squared depth \(8/3\) is in \(K^*/K\): its six coefficient residues satisfy \(a_i+b_i=0\bmod3\), equivalent to belonging to \(\theta K/3K\). In those 729 dual classes take coordinates \(a_i\in\mathbb F_3\); the 252 deepest are exactly \(\sum a_i^2=2\bmod3\). This structural gap, rather than an assumption restricting phases to the dual, is the key ingredient.

## 2. The deepest-label graph and the E8 windows

On any affine dual class, the deepest-edge graph is the additive Cayley graph of \(q(a)=-\sum a_i^2=1\) on \(\mathbb F_3^6\). Its degree is 252 and eigenvalues are \(252,9,-18\) with multiplicities \(1,476,252\), respectively (count the additive characters as in the checker). Its clique size is at most seven: translate a clique vertex to zero. The polar Gram on the other \(k\) vertices is \(I_k+J_k\) over \(\mathbb F_3\). For \(k=7\) this is invertible, contradicting the ambient dimension six. An arbitrary ternary set with deepest differences lies in one affine dual class and has the same cap.

In the usual \(E_8\) coordinates (even-sum integers or half-integers), a nine-point unit root simplex is
\[
0,\quad e_i-e_8\ (1\le i\le6),\quad -e_7-e_8,\quad
(\tfrac12,\tfrac12,\tfrac12,\tfrac12,\tfrac12,\tfrac12,-\tfrac12,-\tfrac12).
\]
All squared pairwise distances are 2. Averaging the count across all translates of this nine-point window gives upper occupancy \(7/9\) whenever every occupied normal root edge must use a deepest horizontal difference and each site has at most one phase. This argument does not depend on periodicity: sum the bound over translates in a large box, and discard just a bounded-width boundary.

For a compact period, \(E_8/2E_8\) has depths \(0,2,4\) in raw units. Indeed the 240 roots use 120 distinct classes (only opposite roots share a class), the polar form of \(q_E(x)=\|x\|^2/2\bmod2\) is nondegenerate by unimodularity, and the Arf sign of the even \(E_8\) gives 120 nonsingular classes. The remaining nonzero classes have depth four from the covering radius 1 of the unscaled \(E_8\) (all can be reached with raw norm at most 4), and parity. The root class graph on 256 vertices is 120-regular, with smallest eigenvalue \(-8\), directly from the \(\pm16\) quadratic character sums in eight variables.

Suppose \(\mu>8/3\) (one phase at each occupied site), \(2<\mu-2s\le8/3\) (every normal root edge becomes deepest), and \(4s<\mu\) (all phases are distinct, since every pair of distinct vertical classes has normal squared depth at most \(4s\)). For \(M\) occupied sites and \(e\) occupied root edges, the normal spectrum yields
\[
2e\ge \tfrac12 M^2-8M.
\]
The labels divide among affine dual classes with counts \(M_b\); each deepest edge stays in a class. Applying the upper nontrivial eigenvalue 9 of the horizontal graph inside each class, and summing, gives
\[
 2e\le\tfrac13\sum_b M_b^2+9M\le \tfrac13M^2+9M .
\]
Here \((252-9)/729=1/3\). Thus \(M\le102\). Importantly, this still applies if different vertical components use different affine dual classes, and if the assignments were never additive.

## 3. Radius and scale, with exact endpoint comparisons

Let \(D(\mu,s)=(\mu/4)^{10}/(27s^4)\). Density equals \(D(\mu,s)\) times average phases per vertical site. In the compact case the average is \(M/256\), and \(s\ge\mu/8\) since each horizontal phase repeats vertically after a distance with squared minimum \(8s\). Always \(\mu\le4\), because an occupied horizontal fiber contains the period \(K\).

* For \(\mu\le2\) any single vertical site is itself a 12D packing with period \(K\). Use center density \(U_{12}<3/50\), a consequence of the [Cohn–Elkies upper bound](https://annals.math.princeton.edu/2003/157-2/p09) with occupied density \(<0.0798\) in [Cohn's table](https://cohn.mit.edu/sphere-packing/), and \(V_{12}=\pi^6/720>(3.14)^6/720\). Its motif occupation is \(\le27U_{12}(4/\mu)^6\). Multiplying by \(D(\mu,s)\), the result is \(\le U_{12}(\mu/(4s))^4\), which is \(<0.124416\) for \(s\ge5/12\). This bound applies at any period.
* For \(2<\mu\le8/3\), each site contains at most seven phases (the gap and clique cap), giving the stated compact fixed-scale cap at \(\mu=8/3,s=5/12\). It is also below the period-free cap for \(s\ge5/8\).
* For \(\mu>8/3\) use uniqueness, and set \(h=\mu-2s\). If \(h>8/3\), no root neighbors survive, so the vertical projection has minimum squared distance \(\ge4s\). The [optimality of \(E_8\)](https://annals.math.princeton.edu/2017/185-3/p07) implies occupancy \(\le1/16\) (the absolute site density is \(\le1/(16s^4)\)). Compact density is at most \(1/27\) by \(s\ge\mu/8\). Period-free with \(s\ge5/8\), use \(D(4,5/8)/16\).
* If \(h\le2\), use full occupancy and \(s\ge(\mu-2)/2\). Compactly this dominates \(\mu/8\) for \(\mu\ge8/3\); the largest endpoint bounds are \(D(8/3,1/3)\) and \(D(4,1)\). Period-free use the maximum of \(5/8\) and \((\mu-2)/2\), switching at \(13/4\). Along the \((\mu-2)/2\) branch logarithmic derivative \(10/\mu-4/(\mu-2)\) has only a minimum, so the two endpoints suffice.
* If \(2<h\le8/3\) and compact, use \(1/27\) if \(4s\ge\mu\); otherwise \(M\le102\) by the spectrum. The scale lower bound \(\max(\mu/8,(\mu-8/3)/2)\) switches at \(32/9\). The density increases to this value, then decreases up to \(\mu=4\), as \(10/\mu-4/(\mu-8/3)<0\) there. This yields the stated compact high-radius cap.
* If \(2<h\le8/3\) and period-free, use occupancy \(7/9\) and \(s\ge\max(5/8,(\mu-8/3)/2)\); the switch is at \(\mu=47/12\), with the same increasing/decreasing check. This gives the stated period-free cap.

All endpoint and comparator comparisons are rational after squaring when necessary. [verify.js](verify.js) replays the entire finite palette histogram, identification of the deep shell, both character tests, the explicit window, and endpoint arithmetic. Run `node verify.js`. The proof rules out the direct reasonable packing candidate and leaves smaller period-free scales and genuinely off-palette constructions open.

## 4. The inherited radius also excludes all positive semidefinite metric stretches

One seemingly natural escape is to lengthen a few normal root directions so that the corresponding pairs no longer need deepest Coxeter–Todd clearance, allowing more occupied sites. This fails for **every** positive-definite normal metric \(G\succeq(2/3)I\), not just scalar or rank-one changes, when the squared minimum is 4. Allow arbitrary assignments and periods as in the first theorem. There is still at most one phase per site. A competitive density would require
\[
\sqrt{\det G}<1/(27\delta_{\mathrm{CHLWW},20}) < (2/3)^{7/2}.
\]
The second strict comparison is checked exactly. Since the seven other eigenvalues are at least \(2/3\), the largest eigenvalue is then less than 1. Consequently **every** root edge still has squared length below 2, so every occupied root edge still requires a deepest difference. The nine-point window therefore gives occupancy at most \(7/9\), and \(\sqrt{\det G}\ge(2/3)^4\) gives density at most \(7/48\), a contradiction. Thus a metric escape at the inherited radius must genuinely shorten some direction below the starting scalar metric, change the vertical support away from a linear image of \(E_8\), or go outside the ternary horizontal alphabet. It cannot be obtained by just stretching the troublesome contacts.
