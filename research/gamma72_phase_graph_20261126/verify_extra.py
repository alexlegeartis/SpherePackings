"""Codegrees, 2x2 rank law, and an independently defined E6-dual candidate."""
import runpy,itertools
from pathlib import Path
import numpy as np
from sympy import Matrix
P=Path(__file__).resolve().parent
d=runpy.run_path(str(P/'verify.py'));R=d['rv'];B=d['flag']
neighbors={int(x):R[B[R^x]] for x in R}
for x,N in neighbors.items():
    counts={}
    for y in N:
        for v in neighbors[int(y)]:
            z=int(v)
            if z!=x:counts[z]=counts.get(z,0)+1
    assert len(R)-1-len(counts)==5713
    assert [sum(v==i for v in counts.values()) for i in (1,2,3,4)]==[94,40,0,2]
    assert max(counts.values())==4
for m,g in enumerate(d['x']):
    p=[(m>>j)&1 for j in range(7)]
    a,b,c,e=p[0]^p[2]^p[3],p[1],p[4]^p[6],p[5]
    assert g==(0 if not m else 4 if not(a or b or c or e) else 5 if (a*e)^(b*c) else 6)
C=2*Matrix.eye(6)
for i,j in ((0,1),(1,2),(2,3),(3,4),(2,5)): C[i,j]=C[j,i]=-1
G=3*C.inv()
assert C.det()==3 and G.det()==243 and all(x.q==1 for x in G)
# z_i^2 <=6*(G^-1)_ii =4.  sqrt(6/3 * 2) = 2 exactly.
H=np.array(G.tolist(),dtype=np.int64);z=np.array(list(itertools.product(range(-2,3),repeat=6)))
M=np.einsum('ij,ij->i',z@H,z);p=(z%2)@(2**np.arange(6))
assert min(M[M>0])==4
roots=set(int(a) for a in p[M==4]);medium=set(int(a) for a in p[M==6])
assert len(roots)==27 and len(medium)==36 and roots.isdisjoint(medium) and roots|medium==set(range(1,64))
links={a:{b for b in roots if a^b in roots} for a in roots}
assert all(len(N)==10 for N in links.values())
assert max(len(links[a]&links[b]) for a,b in itertools.combinations(roots,2))==5
print('PASS maximum projected codegree 4, 2x2 rank law, E6-dual gain 4/3 killed by a required codegree 5')
