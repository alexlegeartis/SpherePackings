# Bacher26 local rigidity and the 27D half-hole contact rank

**No better density found.** The 26D result concerns the original 22+4
period cell and fixed horizontal classes, not all packings of dimension
26. The 27D result concerns a neighborhood of one half-hole over a fixed
horizontal lattice, not all holes or shapes.

The source is [Bacher's original lattice](https://doi.org/10.1007/s002220050180)
in [Martinet's published Gram basis](https://jamartin.perso.math.cnrs.fr/Lattices/bacher.gp).
Use the overlay of the [earlier exact 22+4 decomposition and phase costs](https://github.com/alexlegeartis/SpherePackings/commit/3eab780d02e9f666c1a235c5d1c8578025accb36):
an even rank-26 lattice $L$ of squared minimum 4 and determinant 3, horizontal
rank-22 section $S$ of covolume 8, and normal lattice orthogonally represented
by $D_4=\{n\in\mathbb Z^4:\sum n_i=0\bmod2\}$ with Gram
$D=\operatorname{diag}(6,2,2,2)$. Each label $a$ of the group
$H=\mathbb Z_4^2\times\mathbb Z_2^2$ gives a distinct horizontal residue
of $S$ and a matching normal residue $p_a$ of $D_4$. The squared horizontal
separation of two distinct residues with label difference $k$ is the
known exact value $4-W_k$, with $W_k$ tabulated in the previous note.
The original packing has 64 centers in this product cell, center density
$1/\sqrt3$, and coincides with $L$.

## 26D: independent normal motif motions do not help locally

Vary all of the following continuously: a common horizontal scale
$\sqrt{\lambda}$, a positive-definite normal Gram matrix $Q$, and
an independent displacement $u_a\in\mathbb R^4$ of each of the 64
normal sites. Keep the horizontal residue classes and horizontal lattice
shape otherwise fixed, and keep the normal period $D_4$ as an abstract
period. These parameters do not require the resulting packing to be a
lattice. Modulo a global homothety and a common normal translation, the original data
$(\lambda,Q,u)=(1,D,0)$ give a strict local density maximum in this family.

**Proof.** Normalize squared separation to at least 4, scaling the entire
packing if necessary; that rescaling is part of this family. A literal
norm-4 vector of $S$ is checked by the replay, so $\lambda\ge1$.
Set $u_0=0$ by common translation. Suppose feasible improvements or ties,
distinct from the original after radius normalization and translation, approach
the original configuration. Divide their parameter changes
by their total norm and pass to a unit tangent $(h,E,v)$, where $h\ge0$,
$E$ is a symmetric strain of the normal Gram matrix, and $v_0=0$.
Nonincreasing volume gives
$$
11h+\tfrac12\operatorname{tr}(D^{-1}E)\le0
$$
because the changing part of the cell volume is proportional to
$\lambda^{11}\sqrt{\det Q}$.

For each nonzero label difference $k$, let $z$ be any representative
of the corresponding normal residue with $z^\top D z=W_k$. For every
label $a$, use the two fibers labeled $a$ and $a+k$, choose a horizontal
nearest displacement, and choose the normal period translate realizing $z$.
The infinitesimal separation inequality is
$$
(4-W_k)h+z^\top E z+2 z^\top D(v_{a+k}-v_a)\ge0.
$$
Averaging over $a\in H$ cancels the last term. In particular, the
active vector $e_0/2$ up to sign gives $E_{00}\ge-10h$. Short period
vectors $e_i\pm e_j$ with $1\le i<j\le3$ give
$E_{ii}+E_{jj}\ge2|E_{ij}|$; hence $\sum_{i=1}^3E_{ii}\ge0$.
Therefore the volume derivative is at least $(11-5/6)h$, which
forces $h=0$. The four axial contacts then give $E_{ii}\ge0$,
and the trace inequality makes every diagonal entry zero. The
period vectors force the three small-small entries to vanish.
Mixed active contacts $\frac12(\pm e_0\pm e_i)$ with both sign
products present for each $i=1,2,3$ force the other three entries
to vanish. Thus $E=0$ as well.

The individual contact inequalities are now
$z^\top D(v_{a+k}-v_a)\ge0$ and have zero sum over $a$, so each is an
equality. Taking the finite Fourier transform on $H$ shows, for each
nontrivial character $\chi$ and each $k$ with $\chi(k)\ne1$, that
$z^\top D\widehat v(\chi)=0$ for every active $z$ in that residue.
The checker enumerates the exact 110 active displacements and proves
that the aggregate rows have rank 4 for **every** nontrivial character
(modular rank 4 implies rational rank 4). Thus $\widehat v$ vanishes
away from the trivial character. With $v_0=0$, $v$ is zero.
This contradicts a unit tangent and proves strict local rigidity.

## 27D: no perpendicular shear of the existing Bacher half-hole

The last vector of Martinet's Bacher27 Gram has projected horizontal
class $x=V/2\bmod L$ over the above Bacher26 lattice and squared height 1.
The literal projection equation is checked in the replay. The published
27D lattice minimum 4 implies $d(x,L)^2\ge3$. There are 26 explicit
independent vectors $w$ of squared norm 12 with $w\equiv V\bmod2$,
and hence $d(x,L)^2=3$ and the nearest-point displacements span all
of $\mathbb R^{26}$. Their negatives belong to the same coset.
Consequently the changed adjacent gap at $x+y$ is at most
$$
3+\|y\|^2-\max_w|\langle w,y\rangle|.
$$
Spanning gives a positive lower bound times $\|y\|$ on that maximum.
Thus sufficiently small nonzero $y$ makes the gap strictly smaller
than 3. Even if the shift ceased to be order two, its gain in the
second layer could not help: a first-layer gap smaller than 3 already
requires height larger than 1 to retain squared minimum 4. This
rules out the proposed small fixed-horizontal shear, not larger or
simultaneous deformations.

The historical lattice center-density comparator in both dimensions is
$1/\sqrt3$ ([Nebe–Sloane table](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html)).
The separate 27D unrestricted comparator is [Vardy (1999)](https://doi.org/10.1090/S0002-9947-99-02169-8)
at $1/\sqrt2$; no claim is made against that packing. Ordinary density
is the center density multiplied by the unit-ball volume in the same
dimension.

Exact replay: python3 verify.py, standard library, assertions enabled.
The earlier decomposition supplies the exact horizontal distances; the
replay here independently checks the local metric and phase contact
claims using only integers.
