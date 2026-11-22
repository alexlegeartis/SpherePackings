# Refined 96D lattice bound from 3, 7, and 17 covering centres

The occupied density of a **single 96-dimensional lattice** is at least
**5.64547885516295e-22**. This improves the immediately preceding lattice
bound by **0.2188%**, or **2.6408%** over this repository's pre-dichotomy
lattice bound. The unrestricted nine-translate bound remains
6.25082581067615e-22. This is an existence proof, not an explicit list of
horizontal shears or an optimality claim.

Use the [MW92 order-three dichotomy](../mw92_order3_dichotomy_20261122/):
a deep torus phase either has a sufficiently deep (I−J) multiple, or
forces an order-three phase a₀ of squared depth at least 34/3.
The first branch is unchanged, with its exact 96D Eisenstein Gram.

In the order-three branch use a₀ as the first horizontal shear and the
vertical Gram with diagonal

    (14/3, 5.27266, 5.454682, 5.656144)

and off-diagonal entry half of the diagonal at the earlier of the two
indices. The exact box from the inverse Gram contains every nonzero
vertical difference of squared norm less than 16. Orient each with
positive last nonzero coefficient. This coefficient is always one,
so a short vector at that stage excludes a ball centred at the
negative of the earlier horizontal combination. Since 3a₀ belongs
to MW92, the first coefficient of such a combination is read modulo
three. The complete check finds **3, 7 and 17** distinct excluded
centres at the successive stages after a₀, with required squared
radii at most 10.741, 10.545318 and 10.343856 respectively. Apply the
[Coxeter–Few–Rogers lower bound](https://doi.org/10.1112/S0025579300002059)
to equal-radius coverings at each stage separately. No weighted or
unequal-radius covering hypothesis is used.

The larger of the two vertical Gram determinants is exactly
19889279550909941977/75000000000000000, from the order-three branch.
Multiply it by 47^14 for the squared covolume of the constructed
lattice; its squared minimum is 16, attained in MW92. [Exact results](results.json)
include the comparison with the 96D lattice of
[Bos–Conway–Sloane (1982), Table 1](https://neilsloane.com/doc/Me103.pdf).

Replay from the repository root:

```
python research/mw92_96_equal_cover_20261122/verify.py
```

The arithmetic is rational, including the covering and rounding checks.
The independent `audit_ldl.py` reverse-LDL lattice enumeration checks the same 31 short
vectors up to sign (rather than the inverse-Gram box used here).
