# A non-Nebe 85-dimensional lattice via MW80 specialization

This is an explicit construction with minimum squared norm exactly
16 and covolume 2·41⁸. Its occupied density is at least
3.5478756025421×10⁻¹⁹, about 77.9% of the opening 85D
[repository bound](../README.md). It does not improve either record
and is specified here as a distinct parent for constructive work.

Use the MW lattice L of the curve y²=x³+1+t⁴² over
F₄₁(ω)(t), where ω²+ω+1=0. [Shioda's
construction](https://doi.org/10.2307/2374791) gives rank 80
and covolume 41⁶. The height normalization is the same as the
[characteristic 47 construction](../../mw92_cusp_d4_20261122).
At the 42 simple zeros of 1+t⁴² the fibers are irreducible type II,
and at infinity the scaled equation is y'²=x'³+1+s⁴²,
so the elliptic surface has χ=7 and no reducible fibers.
Every nonzero section has squared height 14+2(P·O), at least 14.
The literal section P=(-1,t²¹) has no poles in any chart and
attains 14.

At t₀=9ω, verified to satisfy t₀⁴²=−1,
specialization into the smooth locus of y²=x³ is the additive
map ρ to F₄₁² with parameter x/y (and 0 at O). As in the
characteristic 47 construction, a finite integral section cannot
pass through the cusp: simultaneous vanishing of x,y would
contradict the simple zero of 1+t⁴². Thus ρ(P)=0 for a nonzero
section forces intersection with the zero section at t₀ and
squared norm at least 16. The two literal sections P and its
image under x↦ωx give the basis c,cω of the additive group,
where c=ρ(P)=9.

Use the vertical Euclidean D₅ with generators
e₀−e₁, e₁−e₂, e₂−e₃, e₃−e₄, e₃+e₄. Define φ:D₅→F₄₁²
by the column pairs (1,0),(26,33),(36,35),(35,37),(17,6)
in the basis c,cω. The map is onto; the standard-library
verifier checks that none of the 2,330 nonzero vertical
vectors of squared length ≤14 have zero image. Thus its
zero-color kernel has squared minimum at least 16.

The packing lattice is the orthogonal-product sublattice
M={(x,z)∈L⊕D₅ : ρ(x)=φ(z)}. Its index over the product
is 41², so its covolume is 41²·41⁶·2. A vector with one
component zero has squared norm at least 16 by the kernel
bounds. A vector with both nonzero has norm at least 14+2.
The matched vector (P,e₀−e₁) attains 16. Hence the center
density is 4^(85/2)/(2·41⁸). The replay uses an outward
Machin interval for π. No sphere decoder is needed.
