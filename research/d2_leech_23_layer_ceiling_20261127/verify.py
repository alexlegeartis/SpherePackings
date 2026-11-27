#!/usr/bin/env python3
"""Exact finite ingredients for the D2 phase ceiling and Golay-covering bound."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import isqrt


def punctured():
    out = [0]
    for shift in range(12):
        out += [x ^ (0xAE3 << shift) for x in out]
    return out


C23 = punctured()
G = [c | ((c.bit_count() % 2) << 23) for c in C23]
assert len(set(G)) == 4096
assert Counter(c.bit_count() for c in G) == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
octads = [c for c in G if c.bit_count() == 8]
five = Counter(t for c in octads for t in combinations([i for i in range(24) if (c >> i) & 1], 5))
assert len(five) == 42504 and set(five.values()) == {1}


def key(p, q):
    return p % 4, q % 4, (p // 4 + q // 4) % 2


C = {}
for p, q in product(range(-4, 5), repeat=2):
    if (p-q) % 2:
        continue
    a = key(p, q)
    n = p*p+q*q
    if a not in C or C[a][0] > n:
        C[a] = n, (p, q)
assert len(C) == 16
assert Counter(n for n, _ in C.values()) == {0: 1, 2: 4, 4: 4, 8: 2, 10: 4, 16: 1}

# Every residue has a representative of raw square <=16, so the above box is complete.
# Points below are projections of actual norm-4 ambient vectors; all sign choices
# on the first two octad positions extend by parity on the other selected positions.
shell = set(product((0, -4, 4), repeat=2))
for octad in octads:
    shell.update(product(*[(-2, 2) if (octad >> i) & 1 else (0,) for i in range(2)]))
for c, j in product(G, range(3)):
    shell.add(tuple((-3 if i == j else 1) * (-1)**((c >> i) & 1) for i in range(2)))
assert all(v in shell for n, v in C.values() if n)


def depth(p, q):
    n = C[key(p, q)][0]
    return 4 - F(n, 8) if n else 0


def anticode(T, d):
    assert len({key(*v) for v in T}) == len(T)
    assert all(depth(p-q, r-s) <= d for (p, r), (q, s) in combinations(T, 2))


anticode([(0, 0)], 0)
anticode([(2, 0), (-2, 0), (0, 2), (0, -2)], 3)
anticode([v for _, v in C.values() if v[0] % 2 == 0], F(7, 2))
anticode([v for _, v in C.values()], F(15, 4))

for i in range(1, 5):
    assert depth(i, i) + F(i*i, 4) >= 4
assert key(4, 4) == key(0, 0)

for size, d, lo, hi in [
    (1, 0, F(7, 2), F(7, 2)),
    (4, 3, F(7, 2), F(37, 10)),
    (8, F(7, 2), F(37, 10), F(379, 100)),
    (16, F(15, 4), F(31, 8), F(4)),
]:
    for mu in (lo, hi):
        assert F(4, size)**2 * (mu/4)**23 <= (mu-d)/4


def low_sqrt(x):
    scale = 10**6
    x = F(x)*scale**2
    return F(isqrt(x.numerator // x.denominator), scale)


lo = F(379, 100)
while lo < F(31, 8):
    hi = min(lo+F(1, 200), F(31, 8))
    a, b, c, e, z = [low_sqrt(lo-d) for d in (F(15, 4), F(7, 2), 3, F(11, 4), 0)]
    assert all((hi/4)**23 <= t*t/4 for t in (b+3*a, c+2*a, e+a, z))
    lo = hi

assert F(23, 24)**23 > F(3, 8)  # The external 23D density is strictly above 1/2.

# Independent finite Golay inputs to the second proof.
Cshort = [c for c in C23 if c.bit_count() % 2 == 0]
assert len(Cshort) == 2**11
assert Counter(c.bit_count() for c in Cshort) == {0: 1, 8: 506, 12: 1288, 16: 253}
assert min(c.bit_count() for c in C23 if c) == 7
assert len(C23) * (1+23+253+1771) == 2**23
heptads = [c for c in C23 if c.bit_count() == 7]
four = Counter(t for c in heptads for t in combinations([i for i in range(23) if (c >> i) & 1], 4))
assert len(four) == 8855 and set(four.values()) == {1}
heavy_octads = [c for c in Cshort if c.bit_count() == 8 and c & 1]
assert len(heavy_octads) == 176
assert all(sum((c >> i) & 1 for c in heavy_octads) == 56 for i in range(1, 23))
assert F(31, 32)**23 < F(1, 2)
print("PASS: Golay design and weighted-cover inputs; exact 16 phase depths and D2 density intervals")
