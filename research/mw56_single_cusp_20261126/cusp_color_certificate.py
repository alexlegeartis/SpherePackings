#!/usr/bin/env python3
"""Five norm-10 color witnesses and one norm-12 kernel witness, all exact."""

from fractions import Fraction as F

P = 29
Q = P * P
T0 = 12


def components(x):
    return x % P, x // P


def encode(a, b):
    return a % P + P * (b % P)


def add(x, y):
    a, b = components(x)
    c, d = components(y)
    return encode(a + c, b + d)


def neg(x):
    a, b = components(x)
    return encode(-a, -b)


def mul(x, y):
    a, b = components(x)
    c, d = components(y)
    return encode(a * c + 2 * b * d, a * d + b * c)


def power(x, n):
    z = 1
    while n:
        if n & 1:
            z = mul(z, x)
        x = mul(x, x)
        n //= 2
    return z


def inv(x):
    assert x
    return power(x, Q - 2)


def trim(f):
    while f and f[-1] == 0:
        f.pop()
    return f


def sub(f, g):
    return trim([
        add(f[i] if i < len(f) else 0, neg(g[i]) if i < len(g) else 0)
        for i in range(max(len(f), len(g)))
    ])


def pmul(f, g):
    r = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            r[i + j] = add(r[i + j], mul(a, b))
    return trim(r)


def peval(f, x):
    z = 0
    for a in reversed(f):
        z = add(a, mul(z, x))
    return z


def prem(f, g):
    f = f.copy()
    u = inv(g[-1])
    while len(f) >= len(g):
        s = mul(f[-1], u)
        j = len(f) - len(g)
        for i, a in enumerate(g):
            f[i + j] = add(f[i + j], neg(mul(s, a)))
        trim(f)
    return f


def pgcd(f, g):
    while g:
        f, g = g, prem(f, g)
    u = inv(f[-1])
    return [mul(x, u) for x in f]


# Coefficients of X=a t^10+b t^5+c and Y=u t^15+d t^10+e t^5+f.
COLORS = [
    ((0, 3, 0, 1, 0, 0, 28), 24, 1),
    ((2, 416, 483, 3, 832, 579, 290), 231, 84),
    ((2, 425, 396, 3, 821, 318, 290), 41, 365),
    ((2, 454, 483, 3, 38, 579, 551), 824, 510),
    ((2, 29, 20, 3, 58, 2, 290), 337, 809),
]


def pullback(params):
    a, b, c, u, d, e, f = params
    return (trim([c] + [0] * 4 + [b] + [0] * 4 + [a]),
            trim([f] + [0] * 4 + [e] + [0] * 4 + [d] + [0] * 4 + [u]))


def check_section(x, y):
    assert len(x) <= 11 and len(y) <= 16
    assert sub(pmul(y, y), pmul(pmul(x, x), x)) == [1] + [0] * 29 + [1]
    a, b = peval(x, T0), peval(y, T0)
    assert a and b
    return mul(a, inv(b))


def main():
    assert pow(2, (P - 1) // 2, P) == P - 1
    assert power(T0, 30) == P - 1 and mul(30 % P, power(T0, 29))
    assert {power(i, 168) for i in range(1, Q)} == {1, 84, 365, 510, 809}
    fifth_powers = {power(i, 5) for i in range(1, Q)}
    covered = set()
    for pars, rho, color in COLORS:
        assert check_section(*pullback(pars)) == rho
        assert power(rho, 168) == color
        covered.update(mul(rho, h) for h in fifth_powers)
    assert covered == set(range(1, Q))

    # Two different polynomial sections; coefficients are low degree first.
    x1 = [1, 0, 60, 0, 289, 0, 262, 0, 85, 0, 28]
    y1 = [812, 0, 403, 0, 182, 0, 657, 0, 652, 0, 377, 0, 248, 0, 466]
    x2, y2 = pullback((362, 502, 337, 834, 307, 213, 36))
    assert check_section(x1, y1) == check_section(x2, y2) == 23
    dx, dy = sub(x1, x2), sub(y1, y2)
    assert len(dx) == 11 and len(dy) == 16  # no meeting at infinity
    assert pgcd(dx, dy) == [17, 1]
    pairing = 5 - 1
    assert 10 + 10 - 2 * pairing == 12

    # Exact normalized endpoint bounds; ordinary and center-density ratios agree.
    r57_cross_sq = F(3, 10) * F(Q * Q, Q * Q - 1) ** 56
    r57_max12_sq = F(3, 2 * Q * Q) * F(6, 5) ** 57
    r58_cross = F(3, 10) * F(Q, Q - 1) ** 28
    r58_max12 = F(3, 2 * Q) * F(6, 5) ** 29
    assert r57_max12_sq < r57_cross_sq < F(548, 1000) ** 2
    assert r58_cross < r58_max12 < F(353, 1000)
    print("PASS: all 840 colors attained at norm 10; kernel norm 12 attained")
    print("PASS: full one-cusp fiber density <0.548 and <0.353 times the 57D/58D incumbents")


if __name__ == "__main__":
    main()
