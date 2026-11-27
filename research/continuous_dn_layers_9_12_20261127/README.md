# Sharp packings by arbitrary real translates of D8–D11 in one extra dimension

**Scope.** With the usual unscaled \(D_m=\{x\in\mathbb Z^m:\sum x_i\equiv0\pmod2\}\), \(m=8,9,10,11\), consider any packing whose centers are a union of full translates of this single fixed-orientation lattice in parallel hyperplanes of \(\mathbb R^{m+1}\). Multiple horizontal translates at the same height are allowed. Horizontal shifts and heights may be arbitrary real numbers and need not be periodic. At squared separation 2 the exact optimal intensities of horizontal layers and their center densities are
\[
I_m=\frac{2}{\sqrt{3-m/4}},\qquad
\delta_m=\frac{2^{-(m+1)/2}}{\sqrt{3-m/4}}.
\]
The values in dimensions 9–12 are respectively \(1/(16\sqrt2),1/(16\sqrt3),1/32,1/32\). The first ties \(\Lambda_9\); none changes an unrestricted density record. In 9D the same conclusion holds for [arbitrary independently rotated full \(D_8\) fibers](ROTATED_D8.md). In dimensions 10–12 the proof also allows two arbitrary horizontal orientations in one packing, but not three or more. The result does not cover partial fibres, and does not claim to solve the continuous cylinder *linear program*.

In particular the 9D result settles the continuous **packing** question for arbitrary phases and independently rotated full \(D_8\) layers left open after Theorem 13 of [Andrew Salmon, Linear Programming Bounds for Fibered Sphere Packings, arXiv:2607.25254v1, §3.2](https://arxiv.org/html/2607.25254v1#S3.SS2), which proves the discrete four-color theorem. The historical comparison for 9D is [Conway–Sloane, What are all the best sphere packings in low dimensions?](https://doi.org/10.1007/BF02574051), and the [9D lattice optimum](https://arxiv.org/abs/2508.20719) concerns a separate class. For 10D use [Best's forty-word code as described by Sloane](https://ems.press/content/book-chapter-files/27252?nt=1), for 12D the [Coxeter–Todd lattice](https://neilsloane.com/doc/Me107.pdf), and for the maintained rows including 11D use the [pinned table](https://github.com/alexlegeartis/SpherePackings/blob/76b46687983161ab92799706f53b85ae3a4ee334/LOWER_BOUNDS_1_96.md).

## 1. Horizontal loss charges the distance to a half-corner

Denote by \(q(x)=d(x,D_m)^2\) the squared horizontal torus distance. Let \(r_i\in[0,1/2]\) be the distance of the \(i\)-th coordinate of \(x\) to \(\mathbb Z\). With
\[
\delta_i=1/2-r_i,\quad d=\min_i\delta_i,\quad
B=\sum_i\delta_i,\quad \epsilon\in\{0,1\},
\]
choose nearest integers in each coordinate and let \(\epsilon\) record their parity. Changing odd parity costs exactly \(1-2\max_i r_i=2d\) by moving one coordinate to its other closest integer. This also covers ties. Consequently
\[
L(x):=m/4-q(x)
   =\sum_i\delta_i(1-\delta_i)-2\epsilon d
   \ \ge\ (1/2-d)B+(m/2-2)d.                 \tag{1}
\]
Here the inequality uses \(\delta_i(1/2-\delta_i)\ge d(1/2-\delta_i)\). Since \(m\ge8\), \(B\ge2\) makes the last bound at least 1. Thus **if \(L<1\), then \(B<2\), and (1) gives \(B\le2L\)**. If \(m>8\), this deduction also holds inclusively at \(L=1\): for \(B>2\) the bound is \(>1\) unless \(d=0\) (when it is \(B/2>1\)) or \(d=1/2\) (when it is \(m/4-1>1\)). For \(B\le2\), equality \(B=2L\) with \(m>8\) forces \(d=0\) in (1) and then \(\delta_i(1/2-\delta_i)=0\) for every coordinate. In this equality case the difference \(x\) is half-integral in every coordinate, with integer coordinates permitted.

The quantity \(B\) is exactly the \(\ell_1\)-distance from \(x\) to its nearest vector with every coordinate in \(\mathbb Z+1/2\).

## 2. Three layers force a normal gap

Let the total squared separation be \(0<\mu\le2\). Consider any three fibers in height order, with successive gaps \(a,b\ge0\) and span \(H=a+b\). If
\[
H^2<3\mu/2-m/4,                               \tag{2}
\]
the three horizontal differences each have
\[
q(x)\ge\mu-H^2>m/4-\mu/2\ge m/4-1;
\]
all three losses are therefore \(<1\). Choose the differences with signs so they sum to zero, and decompose each into an all-half-integral corner plus its shortest coordinatewise error. Three such corners sum to an all-half-integral vector, whose \(\ell_1\)-length is at least \(m/2\). The opposite error sum thus implies
\[
m/2\le B_1+B_2+B_3\le2(L_1+L_2+L_3).
\]
But the packing inequalities give
\[
L_1+L_2+L_3 \le 3m/4-3\mu+a^2+b^2+H^2
             \le 3m/4-3\mu+2H^2 < m/4
\]
by (2), a contradiction.

The cylinder is packed by a locally finite set of phase-height pairs (compact horizontal torus and positive separation). Thus every normal interval of length less than \(\sqrt{3\mu/2-m/4}\) contains at most two fibers, whenever that radical is positive. Iterating the same condition on either subsequence of every other fiber gives the uniform bound \(2T/\sqrt{3\mu/2-m/4}+O(1)\) in every interval of length \(T\). Dividing by horizontal covolume 2 and multiplying by the radius to the power \(m+1\) proves the stated density. The same intensity argument also bounds upper density in Euclidean balls: apply the normal counting bound uniformly in bounded-width bands, use the uniform error for the number of points in each translate of \(D_m\) in a large horizontal ball, and pass to Riemann sums in the normal height after scaling the region. Equivalently one can use normalized large axis-aligned boxes; both recover the displayed density bound.

At other radii let \(t=\mu/2\). The ratio of this analytic upper bound to the displayed sharp value at \(\mu=2\) has square
\[
\frac{(12-m)t^{m+1}}{12t-m}.
\]
It is at most one on \(t\in[17/25,1]\), \([77/100,1]\), and \([87/100,1]\) for \(m=8,9,10\), respectively: the difference \((12-m)t^{m+1}-(12t-m)\) is convex, nonpositive at the stated left endpoint, and zero at 1. For \(m=11\) the estimate is still below the 12D record \(1/27\) for \(19/20\le t\le1\), by the same argument applied to \(729t^{12}-1024(12t-11)\). These ranges are sufficient, not claimed maximal.

## 3. Exact attainment at \(\mu=2\)

Let \(g\) have entries \(1/2\) on the first eight coordinates and zero elsewhere. Let \(x\) have entries \(1/2\) on the first four coordinates and on all coordinates after the first eight, and zero elsewhere. Set \(h=\sqrt{3-m/4}\). At height \(jh\), \(j\in\mathbb Z\), place the two fibers
\[
D_m+jx,\qquad D_m+jx+g.
\]
A half-grid difference with at least one coordinate equal to \(1/2\bmod 1\) has \(q\) equal to one quarter of its number of half-integral coordinates (nearest-integer parity can be corrected on such a coordinate at zero cost). Therefore \(q(g)=2\) and \(q(x)=q(x+g)=(m-4)/4\). The adjacent height term \(h^2\) gives separation exactly 2. For \(m\le10\) all nonadjacent layers have squared vertical separation at least 2. For \(m=11\), exactly lag two needs checking: \(2x\) is an odd integral vector so \(q(2x)=1\), and \(q(2x+g)=2\); the squared vertical term is 1. Larger lags are automatically safe. The centers form a lattice because \(2g\in D_m\) and advancing one height adds the generator \((x,h)\). There are two sites per fiber height and horizontal covolume 2, giving exactly the claimed values. The equality fixture is a verification of attainment, not a priority claim for those lattice configurations.

## 4. What equality says about phase geometry

For a periodic packing attaining the bound at \(\mu=2\), every two-gap span equals \(h\): all spans are at least \(h\), and their average over a height period is \(h\). This assertion also holds for screw periods; only the ordered height list needs to repeat. Therefore the successive gaps alternate between a number \(s\) and \(h-s\). An aperiodic packing with a defect set of zero density need not be tight in every window and is not classified by this argument.

If \(m>8\), the inclusive form of (1) applies to all three losses even at span \(h\). Positive gaps \(s,h-s\) would make \(a^2+b^2+H^2<2h^2\) and again contradict the half-corner sum. Hence the heights are all doubled with gaps \(0,h\). Applying the same equality chain to every consecutive triple gives \(B_i=2L_i\) for each of its differences, so all consecutive horizontal differences are on the half-grid. Their horizontal squared differences are exactly 2 between the coheight pair (half-weight eight), and \((m-4)/4\) across each nearest pair of levels (half-weight \(m-4\)). For \(m=11\) a proposed equality pattern must in addition satisfy the required lag-two squared horizontal depth 1. These conditions, together with the automatically safe longer gaps, are sufficient.

For \(m=8\), if the gaps satisfy \(0<s<1\), every twice-apart phase difference is the odd class in \(\mathbb Z^8/D_8\). To see this, first observe that \(q(p+w)\le1\) whenever \(p\) is integral and \(\|w\|_1\le2\), with equality at \(\|w\|_1<2\) only when \(p+w\) is an odd integer vector. When \(p\) is odd, subtract a signed unit vector in a largest coordinate of \(w\); the squared distance is at most \(1+\|w\|_2^2-2\|w\|_\infty\le1\), strictly smaller if \(0<\|w\|_1<2\). When \(p\) is even, first reduce any magnitude \(>1\) by a signed 2, and if the largest two magnitudes \(A,B\) sum to \(>1\), subtract signed units in both coordinates. After this, every pair of magnitudes sums to at most 1: the transformed pair sums to \(2-A-B<1\), a transformed magnitude \(1-B\) plus any untouched magnitude is at most 1, and the untouched magnitudes have total at most \(2-A-B\). The \(\ell_1\)-norm has not increased. With new largest magnitude \(A\le1/2\) the square sum is at most \(A\|w\|_1\le1\). For \(A\ge1/2\) it is at most \(A^2+(1-A)(2-A)=2A^2-3A+2\le1\). With strict \(\ell_1<2\), equality leaves just one nonzero unit coordinate. All the reductions preserve the parity of the integral translation.

The two consecutive phase differences have \(q>1\) in the positive-gap case, and together differ from an integer vector by an error of \(\ell_1\)-norm at most \(2(s^2+(1-s)^2)<2\). Their twice-apart vertical gap is 1, hence \(q\ge1\) and the equality rule just proved gives the odd class. This means the packing is exactly
\[
D_9\ \cup\bigl(D_9+(w,s)\bigr),\qquad
q(w)\ge2-s^2,\qquad q(w+e_1)\ge2-(1-s)^2
\]
after an overall translation, with \(D_9\) in the usual parity coordinates. These inequalities are sufficient since for each parity of the vertical integer difference the indicated nearest vertical gap is minimal. The four discrete colors are a special case, not assumed here. For instance \(s=1/2\), \(w=(3/8,3/8,1/2,\ldots,1/2)\) gives both horizontal squares \(57/32>7/4\) with genuinely off-palette shifts in the chosen \(D_8\) frame, and a nonlattice equality fixture.

At \(m=8,s=0\), the same-height pair necessarily has squared depth 2, which in (1) forces an all-half-integral shift (there are the two spinor classes). Each doubled level is thus an \(E_8\) layer. The pairwise condition between adjacent integer heights is precisely depth at least 1 in the lattice generated by \(D_8\) and the spinor(s) used on those two levels. This is \(E_8\) if the spinor class is the same, and \(D_8^*\) if the classes differ. All larger height gaps are safe. This condition also characterizes the coincident-height equality patterns without postulating periodic phase order.

Run `node verify.js` for supplemental exact checks of all fixture distances over the full horizontal torus and all relevant vertical lags, the off-palette example, and the rational scale and benchmark comparisons. The upper bound itself is the analytic three-layer argument above.

