# Harada's 51D lattice: a different density tie and exact barriers for arbitrary sheets

26 November 2026. **No density record changed.** The [pinned live table](https://github.com/alexlegeartis/SpherePackings/blob/5af06c1fcb247b511172d52c917306c446de40a1/LOWER_BOUNDS_1_96.md) gives center densities $\delta_n^*=(3/2)^{n/2}/2$ in 51D, 52D and 53D, attained by the [P48n Hurwitz constructions](https://github.com/alexlegeartis/SpherePackings/tree/5af06c1fcb247b511172d52c917306c446de40a1/research/p48n_hurwitz_20261122).

We construct a distinct lattice tie in 51D. A particular hyperplane section has exactly eleven projected phases; their **actual** separation costs yield an exact 51D optimum and a strict 52–53D obstruction within the specified whole-fibre model. Normal sites and assignments may be arbitrary, including nonperiodic. This does not bound noncanonical horizontal offsets or partly occupied horizontal fibres.

## A different 51-dimensional lattice

Take Harada's odd unimodular $H=A_9(C_{51})$ in [arXiv:2611.29304, §4.3 and Theorem 4.4](https://arxiv.org/pdf/2611.29304). We use two published facts: its minimum norm is 5, and the only shortest shadow vectors are $\pm s$, where
$$s=(1^{25},7,-1^{25})/6,\qquad \`s\`^2=11/4.$$
We did not independently enumerate those input facts.

Write a point of $H$ as $w/3$ for an integer vector $w$. The first 25 entries form a word $u$. The next entry satisfies $w_{25}\equiv \sum u_i+3k\pmod 9$. For output coordinate $j$ at index $26+j$, the residue is $\sum_i u_i r_{j-i}+3k$, where $r$ is Harada's 25-entry circulant row and its index is modulo 25. A triangular numerator basis consists of the 25 code rows with identity first block, the row with zero first block and constant 3 thereafter, and the 25 vectors $9e_{26+j}$. Subtracting the first 26 pivots proves that this generates every numerator. The literal triangular determinant is $3\cdot9^{25}=3^{51}$, and division by 3 in each coordinate proves covolume one. Gram products after division by nine are integral; all these identities are rechecked in [verify.py](verify.py).

The even sublattice $H_0$ has index two, hence covolume 2 and minimum at least 6. The first vector in [witness.json](witness.json) attains norm 6. Thus it is an exact **lattice** tie in 51D:
$$\delta(H_0)=(3/2)^{51/2}/2=15463.11943238713\ldots,\qquad \Delta(H_0)=9.34489307985541\ldots\times10^{-10}.$$
This is not the P48n construction under another name: Harada's theta series gives 17,931,672 norm-six vectors, whereas the pinned 51D lattice contains the whole $P_{48n}$, which by itself has [52,416,000 norm-six vectors](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ48a.html). The parent and its theta coefficient are Harada's; no priority is claimed for the even-sublattice operation.

## A full cyclic alphabet of rational true costs

The characteristic vector is $v=2s\in H$, of norm 11. The integer pairing $m=\langle x,v\rangle$ has the parity of the norm of $x$ (check the diagonal congruence on the integral basis). A basis element pairs to one. The vector $v$ is primitive by the parent's minimum, because $11/k^2<5$ for $k\ge2$. The kernel $K=H\cap v^\perp$ is therefore rank 50, even, and has covolume $\sqrt{11}$ by the projected normal spacing $1/\sqrt{11}$.

The horizontal projection of $H$ lands in $K^*$ by integrality. It realizes exactly 11 cosets over $K$, labelled $m\bmod11$: two labels collapse only by multiples of the primitive $v$. Thus it is all of $K^*$. The shadow adds no additional projected colours: $s$ generates the four-element quotient $H_0^*/H_0$ (as $2s\in H\setminus H_0$ and $4s\in H_0$), and $s$ projects to zero.

Reduce $m$ to $[-5,5]$ by multiples of $v$ before taking shortest lengths. The normal squared part is $m^2/11$. For $m=1,3$ the parent norm is at least 5; for $m=0,2,4$ a nonzero point has even norm at least 6. For $m=5$, norm 5 would give $\`s-x\`^2=11/4$, forcing $s-x=\pm s$, contradictory to that pairing. Therefore its odd norm is at least 7. Negative labels follow by negation. The six listed integer numerator witnesses realize these bounds by their exact code membership and Gram costs. This supplies the squared horizontal separations
$$h_0=0,\quad 11(h_1,h_2,h_3,h_4,h_5)=(54,62,46,50,52).$$
The zero cost refers to two copies of an identical phase; within a fibre, nonzero differences have minimum exactly 6.

## All radii and arbitrary normal positions

For $k=1,2,3$ put full translates of $K$ from the eleven phases above individual sites in $\mathbb R^k$. Different phases may occupy the same site if permitted by the distance constraint. The squared packing diameter $\mu$ is at most 6 and the tail separation for phase difference $d$ is at least $\sqrt{\max(0,\mu-h_d)}$. If $I$ bounds the uniform upper number intensity in the tail, the full density ratio is at most
$$\frac{\delta}{\delta^*_{50+k}}\le \frac{2}{\sqrt{11}}(\mu/6)^{(50+k)/2}I.$$
Count complete horizontal period cells in large boxes with uniform vanishing boundary error, then do the same for the normal set. One may pass to large balls by tiling with intermediate-sized cubes. The normal packing bounds hold uniformly by periodically repeating a finite cube with a separating buffer, and the per-colour packing bounds imply local finiteness.

**51D: exact optimum in this model.** For $\mu\le57/10$, one colour has intensity at most $1/\sqrt\mu$, so eleven have ratio at most 0.751175 at the endpoint. For higher $\mu$, every pair of sites is separated by at least $\sqrt{\mu-62/11}$. The square of the resulting ratio is $(4/11)(\mu/6)^{51}/(\mu-62/11)$, whose log derivative has the sign of $50\mu-51(62/11)$. Its maximum is at an endpoint, with value 1 at 6. Equality is realized by spacing $2/\sqrt{11}$, phase $2j$ at site $j$; that is exactly $H_0$, so all separations are inherited. At equality in a finite period every gap has that size and every step is $\pm2$; a sign reversal would repeat a colour after two gaps, less than $\sqrt6$. This identifies the maximizing periodic arrangement up to shifting and reversing, without a uniqueness claim at asymptotic density under sparse defects.

**52D: less than 98% of the live record in this model.** For $\mu\le59/10$, sum the exact hexagonal bound over eleven individual colour sets, obtaining $I\le22/(\sqrt3\mu)$ and ratio at most 0.838617. Above $59/10$, choose a lower endpoint $l\le\mu$ and radius $R=\sqrt{l-52/11}/2$. Only differences $\pm1,\pm2$ can occur within $2R$. Each of those four colours occurs at most once there, because two neighbors of an equal colour would have squared distance less than $4(l-52/11)<l$. Distances to potential neighbors of differences $\pm2$ and $\pm1$ are at least $\sqrt{l-62/11}$ and $\sqrt{l-54/11}$ respectively.

Truncate each Voronoi region to its disk of radius $R$. At most two largest caps of each type are lost. The resulting regions are disjoint in interiors. A cap normalized by $R^2$ has area $2\int_x^1\sqrt{1-t^2}\,dt$ with $x$ the distance to the cutting chord divided by $R$. The verifier certifies an area lower bound on each of 100 rational intervals using the lower endpoint for distance constraints and the upper endpoint for $(\mu/6)^{26}$. Use $\pi>3.14159$ (Machin's formula, alternating series) and a mesh of $10^{-4}$ for the upper left-end cap sums, with every square root rounded up and each $x$ down by exact integer arithmetic. The largest resulting ratio is 0.970969; the certified rational safety threshold is 49/50. No shape, periodicity, colour-balance assumption or discretization of allowed radii remains.

**53D: per-colour Kepler bound.** It gives $I\le11\sqrt2/\mu^{3/2}$, so the ratio is at most $\sqrt{11/27}=0.638285$ at $\mu=6$. Here the tail dimension is three. Sources for the sharp small-dimensional packing inputs: [Cohn's bibliography](https://cohn.mit.edu/sphere-packing/#data) and [Flyspeck](https://doi.org/10.1017/fmp.2017.1).

**Replay:** run `python -B verify.py`. The exact proof verifies the literal integer witnesses and does not rely on the numerical discovery procedure.


[Other normal: 48D and 49D via the 3D fixed space of order 25](FIXED48.md).
