"""Exact continuum and Harper-prefix arithmetic; Python stdlib only."""
from fractions import Fraction as F
from itertools import combinations


def harper5():
    V = sorted(range(32), key=lambda x: (
        x.bit_count(), tuple(i for i in range(5) if (x >> i) & 1)))
    ans = []
    for k in range(3, 17):
        A = set(V[:k])
        N = A | {x ^ (1 << i) for x in A for i in range(5)}
        ans.append(len(N) - k)
    assert ans == [10,11,11,10,12,13,13,12,13,13,12,12,11,10]
    assert 4*(32-min(ans))+192 == 280


def maximum(n, bounds):
    v = [[int(i >= n-k) for i in range(n)] for k in range(n+1)]
    vals = [[f(a) for f in bounds] for a in v]
    cand = []

    def try_(I, c):
        if min(c) < 0:
            return
        assert sum(c) == 1
        a = tuple(sum(ci*v[j][k] for ci,j in zip(c,I)) for k in range(n))
        d = sum((1-t)**2 for t in a)/4 + min(f(a) for f in bounds)
        cand.append((d,a))

    for i in range(n+1):
        try_([i], [F(1)])
    for i,j in combinations(range(n+1),2):
        for u,w in combinations(range(len(bounds)),2):
            d = vals[i][u] - vals[i][w]
            e = vals[j][u] - vals[j][w]
            if d != e:
                r = F(e,e-d)
                try_([i,j], [r,1-r])
    if len(bounds) == 3:
        for i,j,k in combinations(range(n+1),3):
            a,b,c = [vals[t][0]-vals[t][1] for t in (i,j,k)]
            d,e,f = [vals[t][0]-vals[t][2] for t in (i,j,k)]
            det = (b-a)*(f-d)-(c-a)*(e-d)
            if det:
                r = F(-a*(f-d)+d*(c-a),det)
                s = F(-d*(b-a)+a*(e-d),det)
                try_([i,j,k], [1-r-s,r,s])
    opt = max(x[0] for x in cand)
    return opt, sorted(set(a for d,a in cand if d == opt))


harper5()
opt, at = maximum(13, (
    lambda a: a[0]+a[11]+a[12],
    lambda a: sum(a[1:11]),
    lambda a: sum(a[:8])+a[12],
))
assert opt == F(781,196) < 4
assert at == [(F(1,7),)*10+(F(6,7),F(1),F(1))]

opt, at = maximum(12, (
    lambda a: a[10]+a[11],
    lambda a: sum(a[:5]),
))
assert opt == 3 and at == [(F(0),)*12]

print("PASS: Harper prefix boundary, exact all-real cross-gap envelope in 13D and 12D")

