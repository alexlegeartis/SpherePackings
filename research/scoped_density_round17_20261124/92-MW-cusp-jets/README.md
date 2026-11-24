# The MW92 cusp jet filtration, and a six-orbit reduction of a potential 92D improvement

**Scoped results only. No density record changes.** Let (L) be the
Shioda Mordell–Weil lattice of (y^2=x^3+1+t^{48}) over
(𝔽_{47^2}(t)), and set (p=47). Its dimension is 92, determinant
is (47^{14}), and squared minimum is 16;
[Shioda, Theorem 1.1](https://doi.org/10.2307/2374791).
Send a type-II cusp to infinity, so the new equation is
(Y^2=X^3+t^{47}-t). Write (T) for the lattice isometry induced by
(t ↦ t+1), (D=I-T), and (K_j=D^jL).

The **exact squared minima are 18, 20 and 24 at depths 1, 2
and 4**, respectively. At depth 3 the minimum is 22 or 24.
The covolume of (K_j) is (47^{7+2j}). In particular none of these
four lattices or any intermediate lattice between consecutive ones
(K_j,K_{j+1}), (1 ≤ j < 4), improves even the original MW92 density.
For a strictly intermediate subgroup the index acquires a factor of
(p), since the successive index is (p^2).
This does not bound the density of arbitrary sublattices or periodic packings.

## Proof of the jet interpretation and norm bounds

In the local chart (z=1/t), scale by (x=z^{16}X), (y=z^{24}Y).
The equation is (y^2=x^3+z-z^{47}), with simple cusp at (z=0).
Every bad fibre is irreducible; (χ=8) and Shioda's pairing gives
(‖P‖²=16+2(P · O)) for every nonzero section.
Let (V_j) be the subgroup of sections whose local formal parameter
(w=-x/y) vanishes to order at least (j), with the zero section
always included. Then (V_1) is the additive specialization kernel.
The exact coordinate audit checks that specialization is onto
(𝔽_{p^2}^+), and (T) preserves it.
It also checks (Σ(i=0…46) T^i=0) and (T^{47}=I). Thus (T)
has characteristic polynomial (Φ₄₇^2), and
(det D=Φ₄₇(1)^2=p^2). Therefore (K_1=V_1).

Under (T), the parameter becomes ((1+z)^{-8}w(z/(1+z))).
For (w=a z^j+b z^{j+1}+ …), subtracting this series leaves leading
coefficient ((j+8)a) at (z^{j+1}). This also computes the first
nonzero term of the group difference whenever that coefficient is
nonzero: the formal group difference is the ordinary difference
times a unit with constant term one. Two literal polynomial sections
of height 18 give independent (a)'s over (𝔽_p). Iterating
(D) yields independent leading coefficients in each of orders
1 through 39, since (9, …, 46) are nonzero. The coefficient of
(z^j) on (V_j) is a homomorphism with kernel (V_{j+1}), so
each such inclusion has index (p^2). Comparing indices and using
(D V_j ⊆  V_{j+1}) proves (K_j=V_j) for (1 ≤ j ≤ 40).
In fact (V_{40}=V_{41}), because (j+8=0) at (j=39), so
(D V_{39} ⊆  V_{41}). No higher-depth equality or minimum
is claimed. A nonzero vector in (V_j) has squared norm at least
(16+2j), because intersections at any other fibre are nonnegative.
The verifier checks literal integer vectors of norms 20 and 24 in
(D^2L,D^4L). This is important: applying (D) to a height-18 vector
need not attain the bound, because it may introduce other intersections.

## All index-(47^2) sublattices over (D^2 L): only six families remain

The quotient (L/K_2) has order (p^4) and exponent (p): modulo
(D^2), (Σ T^i=pI-C(p,2)D=p(I-23D)=0), and (I-23D)
is invertible there. Choose coordinates ((u,v) ∈ 𝔽_p^2
 × 𝔽_p^2) with (u) the cusp specialization and (D)
acting as ((u,v) ↦ (0,u)). The script constructs this quotient
by a checked linear section of the specialization map and
(pD^{-1}=-Σ(i=0…46) iT^i).
A sublattice (M ⊇  K_2) of index (p^2) defines a two-plane.
Seventeen verified norm-16 vectors and integral isometries of the
parent give witnesses for every nonzero (u). Translations by
powers of (T) fill all (p) possible (v)'s with the same
(det(u,v)).

* If projection to (u) has rank two, the plane is a graph (v=f(u)).
  The value (det(u,f(u))) is a binary quadratic form. All (47^3)
  forms hit at least 68 distinct nonzero (u)'s with a witnessed
  (det(u,v)); finite exhaustive replay included. Every graph
  therefore has squared minimum 16.
* In rank one, the intersection with (v) is a line. If that line
  is not parallel to the projected line, every allowed affine fibre
  crosses the translated line of a norm-16 witness. Again min is 16.
  This exhausts the non-(T)-invariant planes. The remaining planes
  have both lines parallel. For normalizing (u_0=(a,b)), take the
  generators ((u_0,v_0),(0,u_0)); the plane is specified by
  (c=det(u_0,v_0)). Its fibres over (t u_0) have determinant
  (t^2 c), so the same witnesses rule out 1386 of the (48 · 47)
  choices for norm raising. Six more meet literal norm-18 directions.
  The remaining **864 choices are in six explicitly computed
  isometry orbits**, of sizes (12,12,12,276,276,276).
* The rank-zero plane is (K_1/K_2), hence min 18.

Every residual case has min 16, 18, or 20, since (K_2) contains a
norm-20 vector. If any residual representative reaches 20, it gives a
92D lattice denser than the parent by the factor
((20/16)^{46}/47^2 ≈ 12.99). The finite witness subset is
not assumed exhaustive; the six representatives still need a proof
excluding all norm-16 and norm-18 vectors, or counterexamples. Lists of
observed orbit hits by themselves are not a lower-bound certificate.

## Replay

Run `OPENBLAS_NUM_THREADS=1 python verify.py` with NumPy. The frozen
input supplies the exact norm witnesses and the parent and isometry
data. The script independently checks the parent Gram, height pairings
of 92 independent explicit sections, and the literal polynomial action
of (T), with arithmetic in (𝔽_{47}[i]), (i^2=-1).
It verifies subgroup membership by integer divisibility and all graph
and orbit conclusions by finite arithmetic. The `replay.txt` log
lists the residual representatives. For further context see the
[existing MW92 cusp–D4 construction](../../mw92_cusp_d4_20261122/README.md).
