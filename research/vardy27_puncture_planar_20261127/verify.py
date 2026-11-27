"""Exact finite checks for the 27D B*-stock planar flattening proof."""
from collections import Counter
from fractions import Fraction
from itertools import combinations

M = 14
S = 471
LIMIT = (1 << M) - 1


def rot(x, j):
    return ((x << j) | (x >> (M-j))) & LIMIT


def span(basis):
    code = [0]
    for word in basis:
        code += [word ^ x for x in code]
    assert len(code) == len(set(code)) == 2**M
    return code


def check(code):
    hist = Counter(x.bit_count() for x in code)
    assert min(hist.keys() - {0}) == 8
    assert all(i % 2 == 0 for i in hist)
    short = [x for x in code if x.bit_count() == 8]
    assert len(short) == 546
    count = 0
    for i, j, k in combinations(range(2*M), 3):
        patterns = 0
        for w in short:
            f = ((w >> i) & 1) | (((w >> j) & 1) << 1) | (((w >> k) & 1) << 2)
            patterns |= 1 << f
        assert patterns == 255, (i, j, k, patterns)
        count += 1
    assert count == 3276
    # Linear map and eight nonempty fibers imply an equal cardinality.
    assert len(code)//8 == 2048
    return count, sorted(hist.items())


def main():
    cols = [rot(S, j) for j in range(M)]
    rows = [sum(((cols[i] >> j) & 1) << i for i in range(M)) for j in range(M)]
    C = [(1 << j) | (cols[j] << M) for j in range(M)]
    D = [rows[j] | (1 << (M+j)) for j in range(M)]
    assert all((c & d).bit_count() % 2 == 0 for c in C for d in D)
    for label, basis in (("C", C), ("D", D)):
        count, hist = check(span(basis))
        print(label, count, "triples x 8 masks", hist)
    assert 4*49**25 < 3*50**25
    assert 16*3 > Fraction(13, 2)**2
    assert 378-Fraction(25, 2)*32*Fraction(49, 50) < 0
    # Squaring the centre-density ratio coefficient for one unit of total W.
    assert Fraction(2**22 * 8**27 * 2, (2*4**25)**2) == 4
    print("PASS: stock profiles and exact density constants")


if __name__ == "__main__":
    main()
