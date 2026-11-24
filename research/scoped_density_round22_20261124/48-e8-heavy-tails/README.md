# Four exceptional E8 slots: compactifying the heavy-weight domain

In the positive-definite 48-dimensional seven-slot, zero-sum E8
half-period model with weights (a,b,c,d,1,1,1), a<0, and 2≤b≤c≤d,
the following exact tails have center density at most 16,800 for
**all** half-period motifs, linear or nonlinear. This is below the
48D record (3/2)^24 in [Cohn's table](https://cohn.mit.edu/sphere-packing/).
It is not an unrestricted upper bound.

Put D=bcd+a(bc+bd+cd+3bcd). Outside the interval from 1/3500 to
128/1875 the [previous proof](../../scoped_density_round21_20261124/48-e8-four-exceptions-scopes)
already applies. Within each listed determinant interval, any one
weight at or above its threshold suffices; endpoints are inclusive.

| D interval | b threshold | c threshold | d threshold |
| --- | ---: | ---: | ---: |
| 1/3500 to 1/1750 | 12 | 181 | 50,722 |
| 1/1750 to 1/875 | 11 | 128 | 22,510 |
| 1/875 to 2/875 | 9 | 90 | 10,053 |
| 2/875 to 4/875 | 8 | 64 | 4,516 |
| 4/875 to 8/875 | 7 | 46 | 2,033 |
| 8/875 to 16/875 | 6 | 33 | 928 |
| 16/875 to 32/875 | 6 | 24 | 429 |
| 32/875 to 128/1875 | 5 | 13 | 161 |

## Nested scalar periods

Write z=1+3a, F1=b+a(1+3b), F2=bc+a(b+c+3bc).
All three are positive by restricting the positive-definite form.
Consider the integer scalar coefficient system along one common E8
root, with the three unit-weight coefficients equal. With two, one,
or zero of the heavy slots c,d fixed to zero, its rank is 2, 3, or 4,
respectively. The rank-two form has determinant 3F1 and the
rank-three form determinant 3F2, before multiplication by 4; this
follows directly from the rank-one determinant lemma. The period
using just the light slots has squared length 12z.
The exact rank-two and rank-three Hermite constants give available
primitive doubled-period squared lengths with bounds
p₀=12z, p₁²≤64F1, p₂³≤384F2. Shortest scalar vectors are primitive
in their corresponding scalar lattice, and repeating the same
coefficient three times does not affect its gcd.

Each single primitive period of squared length p implies
density ≤2^40(p/8)^24/D^4, by the E8 parity-fibre argument in
the previous proof. In a determinant slab [L,H], choose g to be
the largest thousandth below 4 for which 2^40(g/8)^24≤16800 L^4.
Thus it suffices that any of p₀,p₁,p₂ be ≤g.

## Reducing the unbounded tail to rational intervals

Assume p₀>g and that the stated Hermite upper bounds on p₁ and p₂
are both greater than g. Start from lower bounds b₀,c₀,d₀ that
enforce the particular tail. For s=1/b+1/c+1/d,
z=(s+3D/(bcd))/(3+s), so
1/b+1/c > T=(3g−36H/(b₀c₀d₀))/(12−g)−1/d₀.
Since c≥b and c≥c₀, this bounds b by 2/T and also by
1/(T−1/c₀) when T>1/c₀.

For b in a rational interval [u,v] inside that finite range,
the same inequality bounds c by 1/(T−1/u) when T>1/u.
Also |a|<1/3 implies
F1 = D/(cd)+|a|b(1/c+1/d) ≤ H/(cd)+(b/3)(1/c+1/d).
As F1>g²/64, this yields
c < (H/d₀+v/3)/(g²/64−v/(3d₀))
whenever the denominator is positive. Let C be the smaller
available upper bound. Finally
F2 = (D+|a|bc)/d ≤ (H+vC/3)/d₀, contrary to F2>g³/384
on every remaining interval. The interval width is at most
1/50, and all comparisons in the replay are rational.
Run python3 -B verify.py.
