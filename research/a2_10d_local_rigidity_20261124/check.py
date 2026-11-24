from pathlib import Path
from itertools import combinations

C = sorted(set(tuple(map(int, line.split())) for line in
               Path(__file__).with_name('ternary-18.txt').read_text().splitlines()
               if line.strip()))
assert len(C) == 18 and all(len(c) == 5 and set(c) <= {0, 1, 2} for c in C)
edges = []
for p, q in combinations(range(18), 2):
    D = [i for i in range(5) if C[p][i] != C[q][i]]
    assert len(D) >= 3
    if len(D) == 3:
        edges.append((p, q, D))
for i, j in combinations(range(5), 2):
    assert any(i in D and j in D for p, q, D in edges)
for i in range(5):
    reach = {0}
    while True:
        changed = False
        for p, q, D in edges:
            if i in D and (p in reach or q in reach) and not {p, q} <= reach:
                reach.update((p, q))
                changed = True
        if not changed:
            break
    assert len(reach) == 18
print('PASS', len(edges), 'contact pairs; all ten block pairs; all five graphs connected')
