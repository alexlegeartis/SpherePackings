# 21D: all half-dual phases over the isotropic E7³ period
Date: 2026-11-27. This is a scoped exclusion, not a density improvement.

**Theorem.** Put \(E=E_7\) in its root-length-squared-two convention. Every
periodic packing whose period is \(E^3\) and whose (arbitrarily chosen,
possibly nonlinear) phases belong to \((\tfrac12 E^*/E)^3\) has density
strictly below \(39/50\) of the 21-dimensional packing of Chen, Hu, Li,
Wang and Wu. This includes every attainable radius, and any common scale
of this carrier. It does not include independently scaled factors,
arbitrary real off-palette phases, a different period, or longer-period
label fields within this period.

## Geometry and comparator
Use
\[
E=\{x\in\mathbb Z^8\cup(\mathbb Z+\tfrac12)^8:\ \sum x_i=0\}.
\]
The integer sublattice \(A_7=\mathbb Z^8\cap\mathbf1^\perp\) has
covolume \(\sqrt8\) and index 2 in \(E\), so
\(\mathrm{covol}(E)=\sqrt2\). The minimum squared norm is 2:
integer nonzero sum-zero vectors need at least two unit coordinates,
and half-integer vectors have eight coordinates of absolute value
at least \(1/2\). The vector \(w=(3,3,-1,-1,-1,-1,-1,-1)/4\)
pairs integrally with both \(A_7\) and
\(c=(1,1,1,1,-1,-1,-1,-1)/2\), but \(w\notin E\).
The order of \(E^*/E\) is 2, so \(E^*=E\cup(w+E)\).

Let \(\mu\) denote the minimum squared center separation. Every nonempty
packed motif over \(E^3\) has \(\mu\le2\). For \(M\) distinct phases its
center density is
\[
D(M,\mu)=\frac{M(\sqrt\mu/2)^{21}}{2\sqrt2}
        =\frac M{4096}(\mu/2)^{21/2}.
\]
The [primary comparator, §3](https://arxiv.org/html/2505.02394v1#S3)
is \(D_*=43^{21/2}/(2^{41}3^{23/2})\simeq0.2100363836\).
Center-density comparisons also compare true packing fractions because
both are multiplied by the same 21-ball volume. At \(\mu=2\), volume
alone would allow a gain with \(M\ge861\); mixed half-dual phases enlarge
each factor's alphabet from 128 binary labels to 256 labels, so merely
counting labels does not rule this out.

## Exact 256-phase metric and Fourier table
A phase is represented by \(y/8\), with \(y=4(x+e w)\), \(x\in E\) and
\(e=0,1\). The eight \(y_i\) have a common residue \(3e+2b\pmod4\) for
some \(b=0,1\), and sum zero. Modulo \(E\), keep \(e\) and the \(E/2E\)
coordinates of \(x\). A basis of \(E\) is \(c\) followed by
\(e_i-e_{i+1}\) for \(i=1,\ldots,6\) in zero-based coordinate notation.
Indeed, \(2x_0\) determines the coefficient of \(c\); the remaining
integer vector has coordinate 0 equal to zero and has sum zero, so
it expands in the six differences by successive partial sums. In this
basis \(2w\bmod2E\) is the binary mask \(81\). Accordingly write a
phase as \((a,e)\), \(a\in\mathbb F_2^7\) (encoded as a mask); carries
when both \(e=1\) add mask \(81\).

For each phase, enumerate all \(y\) with \(y^2\le96\) using coordinates
in \([-9,9]\) and the constraints above. The bound on coordinates is
complete since \(10^2>96\). There are exactly 256 residue classes and
each has a representative in this enumeration, so the smallest
\(y^2/64\) observed for each is its true torus squared distance. In
units of \(1/8\), the six depth levels are \(Z=(0,3,4,7,8,12)\) with
sizes \(K=(1,56,63,72,63,1)\). Thus any \(\mu>0\) attained by a motif
is \(t/8\) for an integer \(1\le t\le16\). The minimum across distinct
phases in the product is the sum of their three factor distances,
since their translations can be independently minimized.

All 256 characters are
\[
 \chi_{k,\xi}(a,e)=(-1)^{k\cdot a}\xi^e,
 \qquad \xi^2=(-1)^{k\cdot 81}.
\]
The six distinct rows of character sums over the depth classes are
\[
P=\begin{pmatrix}
1&56&63&72&63&1\\
1&-56&63&-72&63&1\\
1&0&-7&0&7&-1\\
1&-8&-1&8&-1&1\\
1&8&-1&-8&-1&1\\
1&0&9&0&-9&-1
\end{pmatrix}.
\]
They occur \(1,1,72,63,63,56\) times respectively. Importantly, the
replay also checks the reverse averages, not just each row sum:
for every element of depth class \(j\), averaging the characters that
give row \(i\) takes the constant real value \(P_{ij}/K_j\).

## Tensor duals over all radii
For a row triple \(h=(h_1,h_2,h_3)\) average the corresponding product
characters over the six factor permutations. On an element with
depth indices \(s=(s_1,s_2,s_3)\), the value is
\[
F_h(s)=\frac16\sum_{\pi\in S_3}
 \prod_{i=1}^3 P_{h_i,s_{\pi(i)}}/K_{s_{\pi(i)}}.
\]
For every motif \(A\), \(\sum_{x,y\in A} F_h(x-y)\ge0\): this is an
average of squared moduli of character sums. There is no group,
product, linearity, or code-regularity assumption on \(A\).

The seven integer vectors in the companion dual file are all
nonnegative. Let \(R=6\cdot504^3\), so \(RF_h(s)\) is integral. For each
gate \(t_0\) the verifier checks a positive
\[
H=\min_{s:\,Z_{s_1}+Z_{s_2}+Z_{s_3}\ge t_0}
   -R\sum_h Y_h F_h(s).
\]
Each diagonal term contributes \(\sum_h Y_h\). For every other ordered
pair the contribution is at most \(-H/R\); applying the averaged
nonnegativity gives \(M\le1+\lfloor R\sum_hY_h/H\rfloor\).
Checking the 56 sorted triples is exhaustive by permutation symmetry.
The quantitative checks are:

| Actual \(t\) | Bound on \(M\) | Max relative center density in row (approx.) |
| --- | ---: | ---: |
| \(1\)–\(6\) | \(256^3\) | 0.656741 |
| \(7\) | \(256^3/8\) | 0.414234 |
| \(8\)–\(9\) | 274627 | 0.759234 |
| \(10\)–\(11\) | 29906 | 0.679935 |
| \(12\) | 9696 | 0.549644 |
| \(13\) | 4096 | 0.538085 |
| \(14\) | 2645 | 0.756582 |
| \(15\) | 1192 | 0.703591 |
| \(16\) | 663 | 0.770654 |

The \(t=7\) elementary bound is translation averaging of the eight
phases \(0,(e_0-e_j)/2\) for \(1\le j\le7\) in the first factor, zero
in the others. They are distinct mod \(E\), and all pairwise distances
are bounded by squared \(1/2\); every translate meets a \(t=7\) motif at
most once. The lower endpoint of each other row determines the set of
allowed types, while the upper endpoint gives its largest possible
density. Squaring the ratio and multiplying positive denominators
checks rigorously against \((39/50)^2\), without floating-point values.

Replay: run Node on the companion verify.js in this directory. The
finite geometry uses bounded integers, quotient characters have entries
in \(\{\pm1,\pm i\}\), the dual inequalities and comparison use BigInt,
and the printed decimals are informational only. LP optimization is
not part of the certificate.
