"""Exact replay for the all-scale 19D untwisted 2E7 spectral ceiling."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product


def norm4(v):
    """Physical squared norm times four: v holds doubled coordinates."""
    return sum(x * x for x in v)


def label(v):
    assert len(v) == 8 and not sum(v)
    assert len({x % 2 for x in v}) == 1
    coeffs_c = (1, 2, 3, 4, 3, 2, 1)  # doubled A7 prefix coefficients
    k0 = v[0]
    b = k0 & 1
    for j in range(1, 7):
        t = sum(v[:j + 1]) - k0 * coeffs_c[j]
        assert not t % 2
        b |= ((t // 2) & 1) << j
    return b


roots = []
norm_four = []
I = tuple(range(8))
for i, j in product(I, repeat=2):
    if i != j:
        roots.append(tuple(2 * ((k == i) - (k == j)) for k in I))
for U in combinations(I, 4):
    roots.append(tuple(1 if i in U else -1 for i in I))
    for V in combinations(U, 2):
        norm_four.append(tuple(2 if i in V else -2 if i in U else 0 for i in I))
for i in I:
    for V in combinations((j for j in I if j != i), 2):
        v = tuple(3 if j == i else 1 if j in V else -1 for j in I)
        norm_four += [v, tuple(-x for x in v)]

assert len(set(roots)) == len(roots) == 126
assert len(set(norm_four)) == len(norm_four) == 756
assert all(norm4(v) == 8 for v in roots)
assert all(norm4(v) == 16 for v in norm_four)
R = {label(v) for v in roots}
N = {label(v) for v in norm_four}
rad = (3, 3, -1, -1, -1, -1, -1, -1)
assert norm4(rad) == 24
t = label(rad)
assert t != 0 and len(R) == len(N) == 63
assert len(R | N | {0, t}) == 128
assert {r ^ t for r in R} == N

normal_spectrum = Counter(sum((-1) ** ((w & r).bit_count()) for r in R)
                          for w in range(128))
assert normal_spectrum == {63: 1, 9: 28, -1: 63, -7: 36}


def q(v):
    return (2 * (sum(z * z for z in v) + sum(v) ** 2)) % 3


A = list(product(range(3), repeat=6))
assert Counter(map(q, A)) == {0: 225, 1: 252, 2: 252}
assert 225 - 252 == -27  # exact normalized Gauss sign is -1
assert (1 + 6) % 3  # determinant of the nondegenerate polar matrix I_6+J_6
assert (1 + 7) % 3  # a hypothetical 8-clique would have rank-seven polar Gram
S = [v for v in A if q(v) == 1]
horizontal_spectrum = Counter()
for w in A:
    z0 = sum(sum(x * y for x, y in zip(w, v)) % 3 == 0 for v in S)
    k = 3 * z0 - 252
    assert not k % 2
    horizontal_spectrum[k // 2] += 1
assert horizontal_spectrum == {252: 1, 9: 476, -18: 252}

# Spectral transfer inequalities: empty root graph; weighted or injective map.
cap_empty = 7 * 128 // 70
assert cap_empty == 12
assert F(70, 128) - F(1, 3) == F(41, 192)
all_M = range(1, 129)
cap_weighted = max(m for m in all_M if F(41, 192) * m * m <= 16 * m + 126)
cap_injective = max(m for m in all_M if F(41, 192) * m * m <= 16 * m)
assert cap_weighted == 82 and cap_injective == 74


def density_squared(mu, s, M):
    return F(M, 128) ** 2 * (mu / 4) ** 19 / (2 * 27 ** 2 * s ** 7)


mu_c1 = F(8, 3)
mu_c2 = F(32, 9)
assert mu_c1 / 8 == (mu_c1 - 2) / 2
assert mu_c2 / 8 == (mu_c2 - F(8, 3)) / 2
assert mu_c1 < F(19, 6) < 4 < F(38, 9)

squared_bounds = (
    density_squared(F(4), F(1, 2), 12),
    density_squared(F(4), F(2, 3), 82),
    density_squared(mu_c2, mu_c2 / 8, 74),
    density_squared(mu_c1, mu_c1 / 8, 128),
    density_squared(F(4), F(1), 128),
)
ceiling = F(74, 128) * F(8, 27) * F(8, 9) ** 6
assert squared_bounds == (F(1, 36) ** 2, F(82, 128) ** 2 * F(3, 256),
                          ceiling ** 2, F(512, 19683) ** 2, F(1, 1458))
assert max(squared_bounds) == ceiling ** 2

chen_19_squared = F(13 ** 19, 3 ** 18 * 5 ** 21)
assert ceiling ** 2 < F(1, 128) < chen_19_squared
assert density_squared(F(4), F(2, 3), 105) < chen_19_squared
assert density_squared(F(4), F(2, 3), 106) > chen_19_squared

print("PASS: E7 quotient shells, both spectra, occupancy caps and all scale/radius comparisons")
