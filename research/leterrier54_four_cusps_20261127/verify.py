#!/usr/bin/env python3
"""Exact rank-54 four-cusp arithmetic, Python 3 standard library only."""
from fractions import Fraction as F
from math import comb

N, W = 16, 27


def require(b):
    if not b:
        raise ArithmeticError('certificate check failed')


def binom(a, b):
    return comb(a, b) if 0 <= b <= a else 0


def prod(x, y):
    return [sum(x[j] * y[i-j] for j in range(i+1)) for i in range(N)]


one = [1] + [0]*(N-1)
x, y = [0]*N, [0]*N
# The inequalities a^2+ab+b^2 >= (a^2+b^2)/2, and its centered
# analogue, show that [-10,10]^2 covers every term below Q^16.
for a in range(-10, 11):
    for b in range(-10, 11):
        i = 3*(a*a+a*b+b*b)
        j = i + 3*(a+b) + 1
        if i < N:
            x[i] += 1
        if j < N:
            y[j] += 1
require(x[:7] == [1, 0, 0, 6, 0, 0, 0])
require(y[:7] == [0, 3, 0, 0, 3, 0, 0])


def powers(a):
    out = [one]
    for _ in range(W):
        out.append(prod(out[-1], a))
    return out


Xp, Yp = powers(x), powers(y)
P = [prod(Xp[W-j], Yp[j]) for j in range(W+1)]


def h(i, j):
    return sum(binom(W-j,i-k)*2**(i-k)*binom(j,k)*(-1)**k
               for k in range(i+1))


H = [[h(i, j) for j in range(W+1)] for i in range(N)]


def linear(A, v):
    m = len(v)
    T = [[F(i) for i in row] + [F(s)] for row,s in zip(A,v)]
    for i in range(m):
        j = next(j for j in range(i,m) if T[j][i])
        T[i], T[j] = T[j], T[i]
        p = T[i][i]
        T[i] = [a/p for a in T[i]]
        for j in range(m):
            if j != i and T[j][i]:
                p = T[j][i]
                T[j] = [a-p*b for a,b in zip(T[j],T[i])]
    return [r[-1] for r in T]


def jet(origin, C8, n0, n1, n2):
    a = [F(0)]*(W+1)
    for j in range(13):
        target = origin if j == 0 else C8 if j == 12 else 0
        a[j] = F(target-sum(a[i]*P[i][j] for i in range(j)), 3**j)
    for r, nr in enumerate([n0,n1,n2]):
        low = list(range(r, 13, 3))
        high = [j for j in range(13, W+1) if j % 3 == r]
        v = [3*nr*[1,0,0,-6,0][i] - sum(H[i][j]*a[j] for j in low)
             for i in range(5)]
        for j, u in zip(high, linear([[H[i][j] for j in high] for i in range(5)],v)):
            a[j] = u
    f = [sum(a[j]*P[j][i] for j in range(i+1)) for i in range(N)]
    b = [sum(a[j]*H[i][j] for j in range(W+1)) for i in range(N)]
    g = [sum(b[j]*P[j][i] for j in range(i+1)) for i in range(N)]
    require(f[:13] == [origin]+[0]*11+[C8])
    require(g[:5] == [3*(n0+n1+n2),0,0,0,0])
    return [f[13],f[14],f[15],g[5],g[6]]


# Checking the five input basis columns proves all linear formulas.
columns = [jet(*(int(i==j) for j in range(5))) for i in range(5)]
require(columns == [
    [0,0,-398034,5235678,537345900],
    [0,0,-33,27,-567],
    [0,0,15707034,-1194102,24596460],
    [354294,0,0,-91854,268272],
    [0,2499741,0,-367416,4342653],
])
require(columns[0][2] + columns[2][2] == 15309000)
require((columns[0][3]+columns[2][3])//3 == 1347192)


def eval_jet(origin, C8, *ns):
    return [sum(a*b for a,b in zip(row,[origin,C8,*ns])) for row in zip(*columns)]


def gain(s, mu):
    return 3**s * (F(mu)/10)**27


for mu, first in [(8,6),(F(26,3),4),(F(28,3),2),(10,1)]:
    require(gain(first-1,mu)<=1<gain(first,mu))


def counts(rank, minus, dimension):
    # All diagonal nonzero entries can be made +1 except possibly one.
    # In fact this also evaluates any diagonal form with 'minus' entries -1.
    a = [1,0,0]
    for j in range(rank):
        sign = -1 if j<minus else 1
        a = [a[r]+2*a[(r-sign)%3] for r in range(3)]
    return [w*3**(dimension-rank) for w in a]


def witt_restriction(s,e,minus):
    if e%2:
        t = (e-1)//2
    else:
        t = e//2 - int((minus-e//2)%2 != 0)
    return s-e+t


require(counts(1,1,2)==[3,0,6])
require(counts(4,1,4)==[21,30,30])
for ns in ([3,0,6],[21,30,30]):
    require(eval_jet(1,0,*ns)[3] < 0)
require(5*479682 + 4*1660662 == 9*1004562)
require(1004562*729 > 647295138)
require(479682*3**7 > 647295138)
require(eval_jet(1,0,9,0,0)[3] < 0)

# Independent finite case check over all possible restrictions, allowing
# impossible ambient embeddings as well.  Let t count an isotropic
# subspace of the restriction (radical plus its nondegenerate Witt index).
survivors = set()
for s in range(1,26):
    for e in range(s+1):
        for minus in range(e+1):
            ns = counts(e,minus,s)
            t = witt_restriction(s,e,minus)
            require(t>=0)
            for mu in [F(8),F(26,3),F(28,3),F(10)]:
                if gain(s,mu)<=1:
                    continue
                if mu>8:
                    if t>=2 or (mu>F(26,3) and ns[1]) or (mu>F(28,3) and ns[2]):
                        continue
                    if eval_jet(1,0,*ns)[3]<0:
                        continue
                else:
                    gg = eval_jet(1,0,*ns)
                    if gg[4]+21*gg[3]<0:
                        continue
                survivors.add((s,mu,e))
require(survivors == {(1,F(10),0),(6,F(8),0)})
require(eval_jet(1,0,729,0,0)[3:] == [-865264680,18468165240])
require(-F(-865264680,27) == 32046840)
require(F(18468165240,567) == 32571720)
require(F(18468165240-21*865264680,3**7) == 136080)
require(F(136080,21)==6480)
require(F(eval_jet(1,0,3,0,0)[3],9) == 183708)

# A four-shift motif gives the left hand side, while its Fourier
# weights at all vectors of J_10 have the lower bound 1/4.
require(F(eval_jet(1,0,4,0,0)[3],3) == 153090 < F(1347192,4))
require(eval_jet(1,0,5,0,0)[3] < 0)
require(F(5,2)**54 / 3**25 == F(3)*F(5**27,2**27*3**13)**2)
print('PASS: four-cusp interpolation, all rank/type cases, density comparisons and the nonlinear four-phase exclusion.')
