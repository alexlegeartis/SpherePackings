"""Exact finite arithmetic for the 42D exceptional-carrier planar theorem."""

from fractions import Fraction as F


def atan_partial(q, last):
    return sum(((-1) ** j * q ** (2 * j + 1) / (2 * j + 1)
                for j in range(last + 1)), F(0))


def rank_mod3(A):
    A = [row[:] for row in A]
    rank = 0
    for j in range(len(A[0])):
        pivot = next((i for i in range(rank, len(A)) if A[i][j] % 3), None)
        if pivot is None:
            continue
        A[rank], A[pivot] = A[pivot], A[rank]
        pivot_inverse = pow(A[rank][j] % 3, -1, 3)
        A[rank] = [(pivot_inverse * x) % 3 for x in A[rank]]
        for i in range(rank + 1, len(A)):
            old = A[i][j]
            A[i] = [(x - old * y) % 3 for x, y in zip(A[i], A[rank])]
        rank += 1
    return rank


def ij(n):
    return [[1 + (i == j) for j in range(n)] for i in range(n)]


assert rank_mod3(ij(8)) == 7
assert all(sum(row) % 3 == 0 for row in ij(8))
assert rank_mod3(ij(9)) == 9

# Machin's identity, with alternating arctan lower and upper sums.
pi_lo = F(314159, 100000)
pi_hi = F(22, 7)
sqrt3_lo = F(173205, 100000)
assert pi_lo < 16 * atan_partial(F(1, 5), 5) - 4 * F(1, 239)
assert 16 * atan_partial(F(1, 5), 4) - 4 * atan_partial(F(1, 239), 1) < pi_hi
assert 0 < sqrt3_lo and sqrt3_lo ** 2 < 3


t_min, t_clique = F(16, 3), F(14, 3)
e, r_lo, r_hi = F(37, 1000), F(31, 20), F(33, 20)
s_hi = r_hi + e


def a_from_r(r):
    return (t_min - t_clique) / (4 * r ** 2 - 1)


def mu_from_r(r):
    return t_min + a_from_r(r)


def r2_from_mu(mu):
    return (mu - t_clique) / (4 * (mu - t_min))


# Positive cap heights, decreasing h(r,s), site labels distinct, and ring cap.
assert e < r_lo
assert 2 * r_lo ** 2 - 2 * r_lo * e - e ** 2 > 0
assert s_hi < 2 * r_lo
assert 16 * r_hi * e + 4 * e ** 2 < 1
assert 4 * e ** 2 + 4 * s_hi ** 2 * (pi_hi / 11) ** 2 < 1
assert a_from_r(r_lo) < F(2, 25)
assert 4 * s_hi ** 2 * F(2, 25) < t_min


def correction(r):
    return 1 + ((1 + e / r) ** 2 - 1) / 5


def g(mu):
    return (mu / 6) ** 21


def disk_upper(mu, T, M):
    assert mu > T
    return 4 * M * g(mu) / (pi_lo * sqrt3_lo * (mu - T))


def pair(mu):
    return 2 * g(mu) / (3 * (mu - t_min))


mu_L, mu_1, mu_0 = mu_from_r(r_hi), F(10817, 2000), F(540883, 100000)
assert t_min < mu_L < mu_1 < mu_0 < mu_from_r(r_lo) < 6
assert r2_from_mu(mu_L) == r_hi ** 2
assert r2_from_mu(mu_1) < F(1571, 1000) ** 2
assert mu_L > 21 * t_clique / 20  # the clique disk bound is increasing now
base = [
    (F(0), F(15, 4), 0, 6561),
    (F(15, 4), F(41, 10), F(10, 3), 2187),
    (F(41, 10), F(19, 4), 4, 81),
    (F(19, 4), mu_L, t_clique, 9),
]
for lo, hi, T, M in base:
    assert disk_upper(hi, T, M) < 1
    if lo:
        assert disk_upper(lo, T, M) < 1

assert disk_upper(mu_1, t_clique, 9) < correction(r_hi)
assert disk_upper(mu_0, t_clique, 9) < correction(F(1571, 1000))
assert pair(mu_0) < 1
assert pair(F(6)) == 1

# DSW center density squared: lattice determinant 3^7 and squared minimum 6.
assert F(6, 4) ** 42 / 3 ** 7 == F(3 ** 35, 2 ** 42)
print("PASS: phase kernels, annulus geometry, all radius cuts, DSW normalization")
