"""Exact minimum and primitive translation periods for both 11D motifs."""
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
from fractions import Fraction

ROOT = Path(__file__).resolve().parent

@lru_cache(None)
def block_distance(w):
    # Squared distance in doubled coordinates to 4 D4.
    # Any nearest even-sum integer offset lies within one step of the
    # coordinatewise nearest offset: if parity is wrong, change just one.
    nearest = tuple((x + 2) // 4 for x in w)
    return min(sum((x - 4 * h) ** 2 for x, h in zip(w, k))
               for k in product(*(range(x - 1, x + 2) for x in nearest))
               if sum(k) % 2 == 0)

def residue(x, blocks):
    r = [u % 4 for u in x]
    for j in range(0, 4 * blocks, 4):
        r[j] += 4 * (sum(u // 4 for u in x[j:j + 4]) % 2)
    return tuple(r)

def verify(name, blocks, size):
    S = [tuple(map(int, row.split())) for row in (ROOT / name).read_text().splitlines()]
    assert len(S) == size and all(len(x) == 11 for x in S)
    classes = {residue(x, blocks) for x in S}
    assert len(classes) == size
    best = 16
    for a, b in combinations(S, 2):
        w = tuple(x - y for x, y in zip(a, b))
        q = sum(block_distance(w[j:j + 4]) for j in range(0, 4 * blocks, 4))
        q += sum(min(u % 4, -u % 4) ** 2 for u in w[4 * blocks:])
        assert q >= 16, (a, b, q)
        best = min(best, q)
    # Period vectors have minimum 4 in physical squared length.
    assert len(S[0][4 * blocks:]) > 0 and best == 16
    s = S[0]
    shifts = [tuple(y - x for x, y in zip(s, y)) for y in S]
    sym = [h for h in shifts if {residue(tuple(x + t for x, t in zip(c, h)), blocks)
                                  for c in S} == classes]
    assert len(sym) == 1 and sym[0] == (0,) * 11
    V = 2 ** (11 + blocks)
    assert Fraction(size, V) == Fraction(9, 256)
    print(f'{size} sites, min squared 4, primitive covolume {V}, center density 9/256')

if __name__ == '__main__':
    verify('motif144.txt', 1, 144)
    verify('motif288.txt', 2, 288)
