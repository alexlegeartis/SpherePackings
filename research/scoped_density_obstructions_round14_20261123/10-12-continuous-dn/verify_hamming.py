from itertools import combinations
from collections import Counter


def basis(words, n):
    pivots = [0] * n
    for w in words:
        while w:
            i = w.bit_length()-1
            if pivots[i]:
                w ^= pivots[i]
            else:
                pivots[i] = w
                break
    return pivots


def contains(pivots, w):
    while w:
        i = w.bit_length()-1
        if not pivots[i]:
            return False
        w ^= pivots[i]
    return True


n = 10
four = [(i, j, k, l, (1<<i)|(1<<j)|(1<<k)|(1<<l))
        for i, j, k, l in combinations(range(n), 4)]
hist = Counter()
for other_cols in combinations(range(1, 16), n-1):
    cols = (0, *other_cols)
    words = [w for i,j,k,l,w in four
             if cols[i]^cols[j]^cols[k]^cols[l] == 0]
    spans = [basis((w for w in words if w>>i&1), n) for i in range(n)]
    frozen = [w for w in words if all(contains(b,w) for b in spans)]
    H = basis(frozen, n)
    rk = sum(map(bool, H))
    assert rk > 0
    for i in range(n):
        per_coord = basis((w for w in frozen if w>>i&1), n)
        assert sum(map(bool, per_coord)) == rk
    hist[rk] += 1
assert hist == {4:525, 5:4480}
print('all 5005 zero-containing column subsets certified:', dict(sorted(hist.items())))
