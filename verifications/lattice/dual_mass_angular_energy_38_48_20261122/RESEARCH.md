# Approaches for a 1% reduction

Let $U_n$ be the pinned AstroNinja389 upper bound and set the target $T_n=0.99U_n$. At fixed Poisson numerator, a denominator $D$ needs an additional contribution of at least $D/99$. The new dimension 44 certificate adds about $2\cdot10^{-15}D$. Summing the same positive tail values at a few more reciprocal vectors leaves that order-of-magnitude gap. A new Fourier function needs to put appreciable mass at the reciprocal lengths the lattice is forced to have.

## Condition the Fourier sign on a rank-two dual constraint

Set $\Delta(L)>T_n$ and choose a primitive shortest reciprocal vector. In the Fourier squared-length scale let its squared length be $a$ and let $g_2$ be the squared covolume lower bound for any primitive rank-two dual sublattice, obtained from the $n-2$ section bound. For an independent reciprocal vector of squared length $b$, the determinant of the two-vector Gram matrix is at least $g_2$ and at most $ab$. If we put $a$ in a band $I=[a_-,a_+]$, all reciprocal squared lengths then belong to

$$
\bigcup_{k\ge1}k^2 I\ \cup\ [g_2/a_+,\infty).
$$

Fit a separate Gaussian polynomial on each band. Require Fourier nonnegativity only on the allowed intervals and enforce a positive lower value $\eta_I$ on $I$ for the first pair. With the Fourier value at zero normalized to one, optimize the primal numerator against $1+2\eta_I$ and impose the target bound directly. On bounded intervals, an exchange solver can insert the worst continuum sign violation after every solve, then switch to exact Bernstein verification. The maximum over the band bounds and $T_n$ is the new global bound.

The rank-two area lower bounds by themselves give no new forbidden interval above the shortest squared length in ranks 42 and 48. Their existing proofs enclose a reduced reciprocal plane, which is a starting point for conditioning on the whole plane Gram matrix. In a cell where its determinant is at most $B$, apply the $n-3$ section bound to obtain a rank-three squared covolume lower bound $g_3$. Any vector outside the plane has squared length at least $g_3/B$. Below that threshold, enumerate the possible integer combinations in the plane and verify their aggregate contribution, allowing signed Fourier values there. Above it, require a nonnegative tail. The first feasibility calculation is the size of this interval across the Gram cells.

## Use the reciprocal vector to distinguish primal layers

In the unscaled dual choose the same shortest vector $y$. Every primal vector has an integer label $k=\langle x,y\rangle$. Only finitely many labels occur in compact support. The zero layer is a rank-$n-1$ lattice with determinant $\det L\,\|y\|$. If $\|x_i\|^2=a_i$ and the labels are $k_i$, the projected Gram matrix gives the extra inequality

$$
\left(\langle x_1,x_2\rangle-\frac{k_1k_2}{\|y\|^2}\right)^2
\le\left(a_1-\frac{k_1^2}{\|y\|^2}\right)\left(a_2-\frac{k_2^2}{\|y\|^2}\right).
$$

Give each absolute label its own radial feature block inside a single positive-semidefinite harmonic matrix. This preserves harmonic positivity across all vectors and permits separate treatment of the zero layer. On three-line planes the oriented labels satisfy the same signed addition law as the short vectors. The projected Gram inequality and the addition law restrict the off-diagonal sign verification to a finite list of feasible triples for each reciprocal length band. A first model needs only two radial blocks, for the zero layer and its complement.

## Allocate angular energy by plane

The published follow-up keeps the off-diagonal entries of its selected angular subkernel nonnegative individually. For support below squared length two, a rank-two plane has at most three opposite pairs, one per nonzero class modulo twice the plane lattice. The off-diagonal energy can therefore be checked on each plane as a single inequality. Permit a negative pair contribution whenever another pair in the same three-line plane pays for it. This enlarges the feasible set of angular polynomials and keeps the same diagonal gain. The two-line and three-line polynomial domains used by the existing compact-kernel proof provide verifier inputs for this variant. Jointly fit those polynomials with the Fourier coefficients so the spatial and reciprocal gains enter the same target inequality.

Sampled joint Fourier and angular fits produced apparent percent-level gains in several ranks. Exact separation found violations between sampled angles and at the boundary of the three-line domain. The confirmed bounds appear in the [certificate table](README.md). The three proposed extensions target the corresponding sign, geometry, and numerator constraints.
