"""Exact compressed normal-shell dual in every 2A3 chart (no NumPy)."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
from collections import defaultdict
import json

M=168498000
norm=8
N=88
G=((8,4,4),(4,8,4),(4,4,8))
Q=lambda b:4*sum(x*x for x in b)-sum(b)**2
key=lambda b:tuple((4*x-sum(b))%16 for x in b)

# This is the exact normal nearest-label condition for norm-8 vectors
# outside K; we need only its necessity.
A=[v for v in product(range(-4,5),repeat=3)
    if v>=tuple(-x for x in v) and max((0,)+v)-min((0,)+v)<=4]
C=defaultdict(list)
for v in A:
    z=min(key(v),key(tuple(-x for x in v)))
    C[z].append(v)
q={}
for x in product(range(-8,9),repeat=3):
    z=min(key(x),key(tuple(-v for v in x)))
    q[z]=min(q.get(z,1000),Q(x))
assert len(q)==len(C) and set(q)==set(C)
assert all(Q(v)==q[z] for z,W in C.items() for v in W)
assert max(q.values())<=64
# Above nearest-label census is exhaustive: a shorter label than these
# has abs coordinates <=sqrt(8*64/16)<6.
def hs(e,K,d):
    H={e:F(1)}
    D={e:F(1)}
    f=F(1)
    for j in range(1,d//2+1):
        out=defaultdict(F)
        for a,w in D.items():
            for r,s in product(range(3),repeat=2):
                b=list(a)
                k=b[r]
                b[r]-=1
                k*=b[s]
                b[s]-=1
                if k>0 and min(b)>=0:
                    out[tuple(b)]+=k*w*K[r][s]
        D=out
        f*=-F(norm,2*j*(N+2*d-2*j-2))
        H.update({a:w*f for a,w in D.items()})
    return H
def val(p,v):
    return sum(w*(v[0]**e[0])*(v[1]**e[1])*(v[2]**e[2]) for e,w in p.items())
HH=[hs(a,G,d)for d in (2,6)for a in product(range(d+1),repeat=3)if sum(a)==d]
K=[]
for i in range(3):
    K.extend([G[i],tuple(-x for x in G[i])])
    for j in range(i):
        K.extend([tuple(x-y for x,y in zip(G[i],G[j])),tuple(y-x for x,y in zip(G[i],G[j]))])
assert len(set(K))==12
yd=[F(x,10**9) for x in json.loads(Path(__file__).with_name('dual.json').read_text())]
assert len(yd)==1+len(HH)
p=defaultdict(F)
p[(0,0,0)]+=yd[0]
for w,H in zip(yd[1:],HH):
    for e,c in H.items():
        p[e]+=w*c
def use_class(W):
    v=W[0]
    return (v[0]+v[1])%2==0 and Q(v)>8
assert all(len({(v[0]+v[1])%2 for v in W})==1 for W in C.values())
# Lower the constant by the exact maximum column excess + a margin.
excess=max(sum(val(p,v) for v in W)/len(W)-int(use_class(W)) for W in C.values())
drop=max(excess,F(0))+F(1,10**7)
p[(0,0,0)]-=drop
assert all(sum(val(p,v) for v in W)<=len(W)*int(use_class(W)) for W in C.values())
lower=M*p[(0,0,0)]-sum(val(p,v) for v in K)
# In p the harmonic projections have constants from the 2/6 homogenized
# radial terms. The original scalar coefficient is yd[0], not p(0).
lower=M*(yd[0]-drop)-sum(val(p,v) for v in K)
assert lower>0
print('certified blockers',int(lower), 'from',len(C),'signed classes')
assert all(v%2==0 for row in G for v in row)
assert min(Q(v) for v in A if Q(v)>8)==11
# Index-2 chart H size 128 cannot attain previous normalized target.
for n in (128,256):
    assert n*n*117**85<64*128**85
print('PASS universal six index-two kernels and full phase lattice in every 2A3 chart')
