# Exact 41D axial optimum over the exceptional 40D carrier

**No new density record.** Let (L) be the 40-dimensional section in the [actual-distance theorem](README.md). Consider any packing made from full translates of (L) at real heights in one extra coordinate, with every horizontal phase in (L^*/L). The heights may be irregular, several distinct phases may share one height, and the phase sequence need not be additive or periodic. The center density is at most (3^{17}/2^{21}), attained by the nested (E_7)-section construction of [Dorofeev–Sun–Wang, §3.1 and Table 1](https://arxiv.org/html/2607.20359v4).

## Reduce the bound to a one-dimensional interval count

Write (mu) for the squared separation, with (0<mule6) since each horizontal translate contains the known norm-six contacts. If a height interval has length less than (sqrt{mu-T}) for some (0le T<mu), then every pair of layers in it must have actual horizontal squared cross-distance strictly greater than (T). In particular the phases in such an interval are distinct. The exact 40D phase-distance theorem supplies the following caps on their number:

| (T) | cap (M) |
|---:|---:|
| (0) | (6561) |
| (10/3) | (2187) |
| (4) | (81) |
| (14/3) | (9) |
| (16/3) | (1) |

For instance the (81) cap follows from avoiding nonzero isotropic differences, and the (2187) cap from avoiding one cyclic direction whose actual depth is (10/3). Partition the height axis into intervals of length (d<sqrt{mu-T}). Each interval contains at most (M) layers, so upper axial intensity is at most (M/d); let (d) increase to (sqrt{mu-T}). This holds uniformly over starting height and does not require periodicity. Dividing by (operatorname{covol}(L)=81) gives the point intensity in the full space: lattice point counts in horizontal slices have the same leading term for each of the finitely many phases and an error that is uniform in the translate. Multiplying by the 41st power of the sphere radius gives the center-density bound.

The squared ratio of that bound to (3^{17}/2^{21}=(3/2)^{41/2}/(81sqrt{2/3})) is
[
  B_{T,M}(mu)=M^2left(rac{mu}{6}ight)^{41}rac{2/3}{mu-T}.
]
For (T>0) the logarithmic derivative is (41/mu-1/(mu-T)). Its only change of sign is from negative to positive, so the maximum on a closed interval with (mu>T) is at an endpoint. For (T=0) the function is increasing (and tends to zero at the origin).

Take breakpoints (0,15/4,41/10,47/10,107/20,6), using the caps in the table in that order. Exact rational substitution into the three-factor expression gives (B<1) at all positive interval endpoints, except (B_{16/3,1}(6)=1). This proves the upper bound at every radius. The nested 41D section built from (K_7subset K_8) lies in this class: it consists of full translates of (L), and every horizontal projection is in (L^*) by integrality of the ambient lattice. The published determinant (2cdot 3^7) and squared minimum (6) give the density above.

This result excludes an entire real-height search, not arbitrary real horizontal phases or arbitrary 41D packings. The [exact checker](verify_axial41.py) evaluates the ten endpoint comparisons and the reference density normalization.
