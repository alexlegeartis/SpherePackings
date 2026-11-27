"""Rational Delsarte dual + all-radius comparisons for multiple phases per site."""

from collections import defaultdict
from fractions import Fraction as F
from itertools import product

from verify import A, N, R, ceiling, q, t


def ternary_character_sum(v, S):
    # Each set S in this replay is closed under negation.
    n0 = sum(sum(x * y for x, y in zip(v, w)) % 3 == 0 for w in S)
    den2 = 3 * n0 - len(S)
    assert not den2 % 2
    return den2 // 2


def binary_character_sum(v, S):
    return sum((-1) ** ((v & w).bit_count()) for w in S)


def verify_scheme(shells, elements, charsum, expected):
    groups = defaultdict(list)
    for w in elements:
        groups[tuple(charsum(w, sh) for sh in shells)].append(w)
    assert set(groups) == set(expected)
    valencies = expected[0]
    assert tuple(map(len, shells)) == valencies
    for row, frequencies in groups.items():
        for j, shell in enumerate(shells):
            for v in shell:
                assert charsum(v, frequencies) * valencies[j] == row[j] * len(frequencies)


P = ((1, 224, 252, 252), (1, 8, -18, 9), (1, 8, 9, -18), (1, -19, 9, 9))
Q = ((1, 63, 63, 1), (1, -1, -1, 1), (1, 9, -9, -1), (1, -7, 7, -1))
horizontal = [[(0,) * 6]] + [[v for v in A[1:] if q(v) == h] for h in range(3)]
normal = ([0], R, N, [t])
verify_scheme(horizontal, A, ternary_character_sum, P)
verify_scheme(normal, range(128), binary_character_sum, Q)

allowed = {(1, g) for g in (1, 2, 3)} | {(2, g) for g in range(4)} | {(3, 2), (3, 3)}
dual = {(0, 2): 7, (1, 0): F(126, 5), (1, 2): F(882, 5),
        (3, 0): F(112, 5), (3, 2): F(784, 5)}
assert all(v >= 0 for v in dual.values())
for h, g in sorted(allowed):
    total = sum(w * F(-P[i][h] * Q[j][g], P[0][h] * Q[0][g]) for (i, j), w in dual.items())
    assert total >= 1
M = 1 + sum(dual.values())
assert M == F(1944, 5) and int(M) == 388
assert 729 * 18 // (252 + 18) == 48

assert F(388, 128) * F(8, 27) * F(2, 3) ** 6 == F(1552, 19683)
assert F(729, 128) * F(8, 27) * F(1, 2) ** 6 == F(27, 1024)
assert F(729 * 12, 128) * F(8, 27) * F(1, 3) ** 6 == F(1, 36)
assert F(729, 128) * F(8, 27) * F(1, 3) ** 6 == F(1, 432)
upper_squares = [F(1552, 19683) ** 2, F(3, 256) * 7 ** 2 * F(2, 3) ** 12,
                 F(27, 1024) ** 2, F(27, 4096), F(1, 36) ** 2, F(1, 432) ** 2]
assert max(upper_squares) == F(27, 4096) < ceiling ** 2

print("PASS: product character shells, 5 dual weights, and every low-radius comparison")
