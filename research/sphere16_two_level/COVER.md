# Covering the full real torus of an extremal Type II lattice

The **squared covering radius** of any extremal even unimodular lattice
of rank 48 is at most 14/3. In rank 72 it is at most 22/3. In rank
96 the analogous bound 10 is conditional on the existence of an
extremal even unimodular lattice (squared minimum 10). These are
bounds for all real phases, not just finite torsion; exact attainment
of the upper endpoints is not claimed.

Write D=n/12 (respectively 4,6,8). The exact half-period theorem gives
the squared depth of every half class at most D, in half-integer
steps, and a spherical 2-design of contacts at depth D. These facts
are independent of the model of the extremal lattice; see
[the exact half-period proof and arithmetic](https://github.com/alexlegeartis/SpherePackings/blob/6696a479a8a55e440ff8e65f10702702fe7d599c/research/scoped_density_obstructions_round15_20261123/48-72-quarter-torus/HALF_BOUND.md).

Let R denote the squared covering radius and choose a real point x
attaining it, which exists by compactness of the torus. Take a
nearest lattice point w to 2x and put u=2x-w, so that H=||u||²<=R.
Now x=w/2+u/2. Let z run over the relative nearest vectors to the
half-period w/2, of squared depth d0. They occur in antipodal pairs:
if z=w/2-v with v in L, then -z=w/2-(w-v) is also in the same coset.
Consequently

    dist(x,L)^2 <= d0+H/4-max_z |z·u|.

If d0<D, then d0<=D-1/2 by the half-period norm progression, so this
is at most D-1/2+R/4. If d0=D, the 2-design gives mean (z·u)^2 =
D H/n = H/12. Thus the distance bound is D+H/4-sqrt(H/12), a convex
function of H on [0,R]. Its maximum is at H=0 or H=R. In either
case we obtain

    R <= max(D-1/2+R/4, D, D+R/4-sqrt(R/12)).

Put R0=4(D-1/2)/3. In the three stated ranks R0>D and R0>3. If
R>R0 then D-1/2+R/4<R and sqrt(R/12)>1/2, so all three terms in
the maximum are strictly below R. This proves the asserted bounds.

Only the prior half-period arithmetic is needed; the new real-phase
step has no enumeration, approximations, or additional numerical input.

## Finite set of actual Voronoi facet lengths

The same half-class input gives a separate, slightly stronger cutoff on
facet normals than the full covering bound alone. A Voronoi-relevant
vector \(v\) is a shortest vector of its coset \(v+2L\), unique up to
sign. Indeed, if a different vector \(v-2w\) has no larger norm,
the two inequalities for \(w\) and \(v-w\) imply that for \(v\):
the sum of their squared norms is
\((\|v\|^2+\|v-2w\|^2)/2\le\|v\|^2\). Their equality hyperplanes
cannot jointly contain a facet for \(v\), unless they are collinear.
In that case either the same argument gives strict redundancy or
\(w=0,v\), the excluded equal pair. Therefore the half coset at
\(v/2\) has just the two opposite shortest vectors.

At depth \(D\) that set cannot be a spherical 2-design in any of the
stated ranks: a single antipodal pair has rank-one covariance.
Hence every relevant vector has squared norm at most
\(4(D-\tfrac12)\), namely **14 in 48D**, **22 in 72D**, or
**30 in a hypothetical extremal 96D lattice**. This is an upper
cutoff on facet normals, not a classification or an attainment claim.
