#!/usr/bin/env python3
"""Exact replay of the 17D equality witness; standard library only."""
from fractions import Fraction
from functools import cache
from itertools import product
from collections import Counter

COLS = [(1 << (j//4)) | ((j % 4) << 3) for j in range(12)]
SYN = [0] * 4096
for x in range(1, 4096):
    low = x & -x
    SYN[x] = SYN[x ^ low] ^ COLS[low.bit_length()-1]
assert Counter(x.bit_count() for x in range(1,4096) if SYN[x] == 0)[4] > 0
assert all(x.bit_count() >= 4 for x in range(1,4096) if SYN[x] == 0)
assert SYN.count(0) == 128
assert all(SYN[15 << (4*b)] == 0 for b in range(3))

def add(x, y):
    """Packed group addition: bit 5 onwards are half masks."""
    return (x ^ y) ^ SYN[(x >> 5) & (y >> 5)]

def mul(k, x):
    t = 0
    for _ in range(k % 4):
        t = add(t,x)
    return t

@cache
def h4(p):
    a, s = p // 32, p % 32
    return a.bit_count() + 4*min((x & ~a).bit_count() for x in range(4096) if SYN[x] == s)

def chart(symbols, tag):
    a = sum((15 << (4*b)) for b in range(3) if symbols[b] & 1)
    s = sum(1 << b for b in range(3) if symbols[b] & 2)
    return 32*a + s + 8*tag

for syms in product(range(4), repeat=3):
    for tag in range(4):
        expected = 4*sum(t != 0 for t in syms) if any(syms) else 8*(tag != 0)
        assert h4(chart(syms,tag)) == expected

def core(z):
    assert sum(z) % 2 == 0
    symbol = 0
    for j in range(4):
        if z[j] % 2:
            symbol ^= j
    tag = (sum(z)//2 % 2) ^ (2*(z[3] % 2))
    return chart([symbol]*3,tag)

ETA = 26497

def phase(y):
    a,b,c,d,e = y
    return add(core((a,b,c,d-e)), mul(e,ETA))

counts = Counter()
ties = Counter()
for y in product(range(-2,3), repeat=5):
    if sum(y) % 2 or not any(y):
        continue
    q4 = 2*sum(c*c for c in y[:4]) + 4*y[4]**2
    if q4 >= 16:
        continue
    assert q4 + h4(phase(y)) >= 16, (y, q4, h4(phase(y)))
    counts[q4] += 1
    if q4 + h4(phase(y)) == 16:
        ties[q4] += 1
assert sum(counts.values()) == 320
assert ties[6] > 0
print("rank-five Hamming parity check, size 128; chart metric checked for all 256 phases")
print("short nonzero normal vectors by 4 x norm squared:", dict(sorted(counts.items())))
print("short normal vectors attaining equality:", dict(sorted(ties.items())))
print("minimum squared = 4 (including H period 2e_i), center density =", Fraction(1,32)/Fraction(1,2))
