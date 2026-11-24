# Five cusp polynomials rule out the proposed 96D second-jet lift

This is not a new density bound. It closes the explicitly proposed
cusp-derivative construction for every choice of its one-dimensional
discarded derivative subspace.

Use the MW92 equation and height normalization in the earlier
[cusp + D4 construction](../mw92_cusp_d4_20261122/README.md),
with the type-II cusp at t0=18+2i over F47[i], i²=−1. In the change
of variable t=t0(v+i)/(v−i), first multiply the x,y coordinates by
(v−i)^16,(v−i)^24, respectively, and then divide by k²,k³,
where k=9+9i and k⁶=−2i. The normalized equation is
```
Y²−X³ = v^47−v.
```
Here (v−i)^48−(v+i)^48=−2i(v^47−v). An exact solution with
degrees 18 and 27 yields a section of squared height 18:
it has exactly one simple meeting of the zero section, at
t0, as seen from the degree excess (2,3) there. It is integral
in every other finite chart and at the original infinity
(use its minimal coordinate changes of degrees 16 and 24).
Shioda's height is 16+2(P·O) on this curve. See
[the original height-pairing paper](https://doi.org/10.14992/00009994).

The slope of the derivative of x/y at t0 is read from
`(X18/Y27)/(k·2i·t0)`. The additive formal parameter has the
same derivative at the zero section: its first-order change of
coordinate is proportional to (t−t0)(x/y)^7, which vanishes to
higher order there.

The five pairs in `witness_polys.json` satisfy the displayed
equation coefficientwise. In the original cusp basis, their
literal projective derivative lines are 0, 1, 2, 3 and 16.
The last two previously unresolved orbits are represented by
0 and 16. For dilation v→5v, take λ=37+40i with λ⁶=5 and scale
the two outputs by λ⁻²,λ⁻³, respectively. This multiplies the
leading ratio by λ·5⁻⁹=35+i. Together with coefficient conjugation
it covers all 48 projective derivative lines from these five
polynomials. Run `python -B verify.py` to check the exact field
and polynomial arithmetic, derivative values, and disjoint orbit
coverage. No lattice search is needed to replay the obstruction.

For the proposed 96D graph, the horizontal component is required
to have zero cusp residue, and its derivative is tested by one
nonzero F47-linear functional. That functional has a projective
kernel line, and the corresponding section above is a pure
horizontal vector of squared norm 18. Thus the hoped-for squared
minimum 20 fails. With the proposed tail
`sqrt(10)(D4 + Z (48,2,4,10)/47)`, the actual minimum is 18:
all horizontal norms in the cusp kernel are at least 18, the
pure vertical colour-zero kernel has minimum squared 20, and
every mixed vector has squared norm >20. This only rules out
this derivative functional and tail strategy, not other 96D
constructions.
