"""Exact integer shell and complement checks for the rank-five packet reduction."""
from collections import Counter
from itertools import product
from math import isqrt
import numpy as np

e=((0,1),(1,2),(2,3),(3,4),(2,5))
G=np.eye(6,dtype=np.int64)*2
for a,b in e:G[a,b]=G[b,a]=-1
Ginv3=np.asarray([[4,5,6,4,2,3],[5,10,12,8,4,6],[6,12,18,12,6,9],
                   [4,8,12,10,5,6],[2,4,6,5,4,3],[3,6,9,6,3,6]],dtype=np.int64)
assert np.array_equal(G@Ginv3,3*np.eye(6,dtype=np.int64))
B=[isqrt(4*int(x)) for x in Ginv3.diagonal()]
Z=np.asarray(list(product(*(range(-b,b+1) for b in B))),dtype=np.int64)
L=np.asarray(list(product(range(3),repeat=6)),dtype=np.int64)
norm=np.sum((Z@G)*Z,axis=1)
p=3**np.arange(5,-1,-1)
d=np.full(729,100,dtype=np.int64)
np.minimum.at(d,(Z%3)@p,norm)
assert Counter(d)=={0:1,2:72,4:270,6:240,8:144,12:2}
r=np.array([2,1,0,2,1,0],dtype=np.int64)
assert np.all((r@G)%3==0)
W=L[:243]
ty=np.sum((W@G)*W,axis=1)//2%3
ty[1:][ty[1:]==0]=3
lift=(W[:,None,:]+np.arange(3)[None,:,None]*r)%3
ds=d[lift@p]
assert {t:{tuple(sorted(int(i) for i in row)) for row in ds[ty==t]} for t in range(4)}=={0:{(0,12,12)},1:{(2,8,8)},2:{(4,4,4)},3:{(6,6,6)}}
# Only these quotient-type pairs have exactly one forbidden radical offset;
# all other nonempty forbidden sets for distinct labels are affine lines.
BADS={(0,1),(1,0),(1,1)}
GOOD={(1,3),(2,2),(2,3),(3,1),(3,2),(3,3)}
LINE={(0,2),(0,3),(1,2),(2,0),(2,1),(3,0)}
grid=np.asarray(list(product(range(3),repeat=2)))
for t,u in product(range(4),repeat=2):
    if t==u==0:continue
    for v in np.unique(ds[ty==t],axis=0):
        for w in np.unique(ds[ty==u],axis=0):
            F=grid[(v[:,None]+w[None,:]).reshape(-1)<8]
            if (t,u) in BADS: assert len(F)==1
            elif (t,u) in LINE: assert len(F)==3 and not (F.sum(0)%3).any()
            else: assert (t,u) in GOOD and len(F)==0
# Fix the canonical depth-0 or depth-2 minimal lift of every type 0/1 label.
canonical=lift[np.arange(243),ds.argmin(1)]
wr= (W[:,None,:]-W[None,:,:])%3
idx=wr@p
q=np.flatnonzero(ty<=1)
ij,jk=np.where(ty[idx[np.ix_(q,q)]]<=1)
i,j=q[ij],q[jk]
assert np.array_equal((canonical[i]-canonical[j])%3,canonical[idx[i,j]])
# Every shape of size >=4 has >=2 projected elements on every one-dimensional quotient.
for b in range(512):
    Q=grid[[(b>>u)&1==1 for u in range(9)]]
    if len(Q)>=4:
        for normal in [(1,0),(0,1),(1,1),(1,2)]:
            assert len(set(Q@normal%3))>=2
print("PASS: all phase-depth masks, canonical affine offsets, and noncollinearity checks")
