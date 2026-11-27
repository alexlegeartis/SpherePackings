#!/usr/bin/env python3
"""Exact truncated q-series and rational interpolation; standard library only."""
from fractions import Fraction as Q
from math import comb

N = 4


def pad(xs):
    return list(xs) + [0] * (N + 1 - len(xs))


def mul(a, b):
    return [sum(a[j] * b[i-j] for j in range(i+1)) for i in range(N+1)]


def pw(a, k):
    out = pad([1])
    for _ in range(k):
        out = mul(out, a)
    return out


def eta_product(modulus, exponent):
    out = pad([1])
    for d in range(modulus, N+1, modulus):
        p = pad([])
        for j in range(N//d+1):
            if exponent >= 0:
                p[j*d] = (-1)**j * comb(exponent, j) if j <= exponent else 0
            else:
                p[j*d] = comb(-exponent+j-1, j)
        out = mul(out, p)
    return out


def solve(A, b):
    m = [list(map(Q, row + [v])) for row, v in zip(A, b)]
    for j in range(len(m)):
        i = next(i for i in range(j, len(m)) if m[i][j])
        m[i], m[j] = m[j], m[i]
        p = m[j][j]
        m[j] = [x/p for x in m[j]]
        for i in range(len(m)):
            if i != j:
                p = m[i][j]
                m[i] = [a-p*c for a, c in zip(m[i], m[j])]
    return [row[-1] for row in m]


T = pad([])
for a in range(-4, 5):
    for b in range(-4, 5):
        n = a*a+a*b+b*b
        if n <= N:
            T[n] += 1
# The omitted region has norm > 4, since n >= 3 max(a²,b²)/4.
C = mul(eta_product(1, 9), eta_product(3, -3))
D = mul(pad([0, 27]), mul(eta_product(3, 9), eta_product(1, -3)))
assert pw(T, 3) == [a+b for a, b in zip(C, D)]

B = [mul(pw(T, 20-3*j), pw(D, j)) for j in range(7)]
W = [mul(pw(T, 20-3*j), pw(C, j)) for j in range(7)]
rows = [[p[i] for p in B] for i in range(3)] + [[p[i] for p in W] for i in range(4)]
values = []
for r in (0, 1):
    coeff = solve(rows, [1, 0, 0, r, 0, 0, 0])
    values.append((sum(c*p[3] for c, p in zip(coeff, B)), sum(c*p[4] for c, p in zip(coeff, W))))
assert values == [(-1500, 360855), (0, 360360)]

anti = [mul(pw(T, 20-3*j), pw([d-c for c, d in zip(C, D)], j))[1:4] for j in (1, 3, 5)]
assert anti == [[-66, -1170, 35034], [42, 1422, -32034], [150, -7650, 110850]]
a, b, c = anti
det = sum(a[i]*(b[(i+1)%3]*c[(i+2)%3]-b[(i+2)%3]*c[(i+1)%3]) for i in range(3))
assert det == -1889568000

assert 360855 == 495 * 3**6 and 1500*(3**6-1) == 1092000
for k in (1, 2):
    assert Q(22, 24)**20 < Q(3**(26-k), 2**40) < 1
assert Q(2**20, 3**16) < 1 < Q(2**40, 3**25)
print("PASS: exact jets, [q³]Theta=1500(r−1), [q⁴]Fricke=495(729−r), density comparisons")
