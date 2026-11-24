"""Exact 12D nonlinear ternary exclusion at two sheared cone masks."""
from fractions import Fraction as F
from itertools import product
from math import isqrt,lcm
import subprocess
from functools import lru_cache
import numpy as np

# E6 root Gram with simple-root labels 0--1--2--3--4, and 2--5.
G=2*np.eye(6,dtype=np.int64)
for i,j in [(0,1),(1,2),(2,3),(3,4),(2,5)]:
    G[i,j]=G[j,i]=-1
I=np.array([[4,5,6,4,2,3],[5,10,12,8,4,6],[6,12,18,12,6,9],
    [4,8,12,10,5,6],[2,4,6,5,4,3],[3,6,9,6,3,6]],dtype=np.int64)
assert np.array_equal(G@I,3*np.eye(6,dtype=np.int64))
invd=[F(int(v),3) for v in np.diag(I)]
N=3**12; pw=3**np.arange(12)
V=(np.arange(N)[:,None]//pw)%3
# Orbit partition under diagonal Weyl group and joint negation.
maps=[]
for i in range(6):
    U=V.copy()
    U[:,i]-=U[:,:6]@G[:,i]
    U[:,i+6]-=U[:,6:]@G[:,i]
    maps.append((U%3)@pw)
maps.append((-V%3)@pw)
done=np.zeros(N,bool); O=[]
for r in range(N):
    if done[r]: continue
    Q=[r];done[r]=True
    for s in Q:
        for m in maps:
            z=int(m[s])
            if not done[z]:done[z]=True;Q.append(z)
    O.append(Q)
assert len(O)==124 and O[0]==[0]

@lru_cache(None)
def shell(r,H):
    L=[isqrt(int(H*u)) for u in invd]
    S=[[j for j in range(-h,h+1) if j%3==s] for s,h in zip(r,L)]
    U=np.array(list(product(*S)),dtype=np.int64).reshape(-1,6)
    Q=np.sum((U@G)*U,axis=1)
    return U[Q<=H],Q[Q<=H]

for a,c,t,nnb,nrep in [(F(3,2),F(3,2),F(3,2),5654,14),
    (F(36,25),F(81,50),F(36,25),5940,15)]:
    den=lcm(a.denominator,c.denominator,t.denominator)
    # Completing one square gives the Schur lower bound on the other
    # slot. Then Cauchy in G gives the coordinate-wise inverse-Gram bounds.
    H1=int(18/(a-t*t/(4*c)))
    H2=int(18/(c-t*t/(4*a)))
    adj=[]; tight=set()
    for Q in O[1:]:
        v=V[Q[0]]
        X,x=shell(tuple(v[:6]),H1)
        Y,y=shell(tuple(v[6:]),H2)
        Z=int(a*den)*x[:,None]+int(c*den)*y[None,:]+int(t*den)*(X@G@Y.T)
        if np.any(Z<18*den):continue
        adj.append(Q)
        if np.any(Z==18*den):tight.update(Q)
    neigh=sorted(u for Q in adj for u in Q)
    assert (len(neigh),len(adj))==(nnb,nrep)
    S={
        F(3,2): [0,140,102200,102256,142884,143080,370297,375416,
            387826,395772,410366,416044,451084,477364,487590,436484,
            430641,472170,506573,518160,530573],
        F(36,25): [0,3787,102200,103707,137786,143080,365000,465010,
            499320,198753,222753,239384,256913,297784,339313,375416,
            398686,416044,435030,475904,518160]
    }[a]
    diffs=set(map(int,(((V[S][:,None]-V[S][None,:])%3)@pw).ravel()))
    assert len(set(S))==21 and diffs-{0}<=set(neigh) and diffs&tight
    # With both integral slots nonzero, Cauchy bounds norm by
    # 2(a+c-t); a pure slot has norm at least 2 min(a,c).
    assert min(a,c,a+c-t)>=1
    target=22
    print('Sheared ternary mask',a,c,t,'target',target,flush=True)
    raw=' '.join(map(str,[target,nnb,nrep]+neigh+[Q[0] for Q in adj]))
    subprocess.run(['./cliques'],input=raw,text=True,check=True)
print('Both sheared nonlinear ternary motif capacities are exactly 21.')
