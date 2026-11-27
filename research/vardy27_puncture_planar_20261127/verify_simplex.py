"""Replay every puncture witness for the per-codeword phase bound."""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
from verify import M, S, rot, span


def main():
    cols = [rot(S, j) for j in range(M)]
    rows = [sum(((cols[i] >> j) & 1) << i for i in range(M)) for j in range(M)]
    basis = [[(1 << j) | (cols[j] << M) for j in range(M)],
             [rows[j] | (1 << (M+j)) for j in range(M)]]
    codes = [set(span(b)) for b in basis]
    triples = set(combinations(range(2*M), 3))
    covered = [set(), set()]
    sizes = [0, 0]
    path = Path(__file__).with_name("simplex_witnesses.txt")
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        b, x, y, z = map(int, line.split())
        assert b in (0, 1)
        nonzero = {x, y, z, x^y, x^z, y^z, x^y^z}
        assert len(nonzero) == 7
        assert all(w.bit_count() == 8 and w in codes[b] for w in nonzero)
        sizes[b] += 1
        columns = [(x >> i & 1) | ((y >> i & 1) << 1) | ((z >> i & 1) << 2)
                   for i in range(2*M)]
        for i, j, k in triples - covered[b]:
            a, d, e = columns[i], columns[j], columns[k]
            if a != 0 and d not in (0, a) and e not in (0, a, d, a^d):
                covered[b].add((i, j, k))
    assert all(c == triples for c in covered)
    assert len(triples) == 3276 and sizes == [214, 241]
    assert 48*Q(7, 8)**2 > 36
    assert Q(147, 256) > Q(39, 40)**25
    assert 378-Q(25, 2)*32*Q(39, 40) < 0
    assert Q(7, 8)*(7+Q(10, 8)) > 7
    print("PASS:", sizes, "literal constant-weight subspaces cover all",
          len(triples), "coordinate triples in each stock")


if __name__ == "__main__":
    main()
