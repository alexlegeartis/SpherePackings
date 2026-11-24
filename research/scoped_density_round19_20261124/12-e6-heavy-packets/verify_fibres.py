"""Actual single-coordinate fibres in the equal-scale 12D mask."""
from itertools import combinations
from collections import Counter
from pathlib import Path
import sys,runpy

path=(Path(sys.argv[1]) if len(sys.argv)>1 else
      Path(__file__).resolve().parents[2]/
      "scoped_density_round18_20261124/12-e6-ternary/replay.py")
r=runpy.run_path(str(path))
d,I,V=r["d"],r["index"],r["V"]
assert V[0]==(0,)*6
Ch=[i for i in range(1,len(V)) if d[i]>=8]
assert len(Ch)==146
def diff(i,j):
    return I[tuple((x-y)%3 for x,y in zip(V[i],V[j]))]
E={i:{j for j in Ch if j!=i and d[diff(i,j)]>=8}for i in Ch}
out=Counter()
for i,j in combinations(Ch,2):
    if j not in E[i]:
        continue
    assert not (E[i]&E[j])  # no fourth site after translating one to zero
    shape=tuple(sorted((d[i],d[j],d[diff(i,j)])))
    col=V[j]==tuple(-x%3 for x in V[i])
    out[col,shape]+=1
assert out=={(False,(8,8,8)):1440,(False,(8,8,12)):216,
            (True,(8,8,8)):72,(True,(12,12,12)):1}
print('PASS: exact actual-coordinate fibre size three, with anchored triple types',out)
