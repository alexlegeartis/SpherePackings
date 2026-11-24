"""Literal integer audit. The analytic identification of rho is separately pinned."""
import gzip,json
from itertools import product
from math import prod
from pathlib import Path
def read(path):
 with gzip.open(Path(__file__).with_name(path),'rt')as f:return json.load(f)
G=read('gram.json.gz')
R=read('colour_q16_witnesses.json.gz')
assert len(G)==len(R['characters'])==92 and R['p']==47 and R['t']==112
T=list(zip(*G)); J=list(zip(*R['characters']))
assert len(J)==2 and 0 not in map(len,T)
M=set()
for v in R['vectors']:
 assert len(v)==92 and all(type(z)is int for z in v)
 y=[sum(a*b for a,b in zip(v,c))for c in T]
 assert sum(a*b for a,b in zip(v,y))==16
 lab=tuple(sum(a*b for a,b in zip(v,c))%47 for c in J)
 assert lab not in M;M.add(lab)
assert M==set(product(range(47),repeat=2))-{(0,0)}
print('literal q16 witnesses cover exactly',len(M),'nonzero values at cusp 18+2i')
