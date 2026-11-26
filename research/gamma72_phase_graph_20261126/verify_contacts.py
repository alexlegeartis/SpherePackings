"""Exact marked shell difference, hence old and new full lattices are not isometric."""
import runpy,json,itertools
from pathlib import Path
from functools import lru_cache
from math import isqrt
import numpy as np
from sympy import Matrix
P=Path(__file__).resolve().parent
d=runpy.run_path(str(P/'verify.py')); M=d['M']; pc=d['pc']; par=d['d']; B0=d['other']; raw=d['raw']
av=np.arange(4096,dtype=np.int32)
S=np.zeros((4096,4096),dtype=np.uint8)
for i,aa in enumerate(M):
    if aa:S[i,list(aa)]=1
lists=[np.array(sorted(x),dtype=np.uint16) for x in M]
def shell(a,b,l):
    if l==4:return 2*S[a,b].astype(np.int32)
    if l==6:return 2*pc[a&b].astype(np.int32)
    if l==8:return 48*((((a!=0)|(b!=0))&(S[a,b]==0)&(pc[a&b]==0)).astype(np.int32))
    return 0
@lru_cache(None)
def count(a,b,t,h):
    need=2*h;n=0;rel=(a,b,0)
    for i in range(3):
        j,k=rel[(i+1)%3],rel[(i+2)%3]
        if j==k:n+=int(shell(rel[i]^j,np.int32(t),need))
        u,v=rel[i]^j,rel[i]^k
        n+=sum(int(np.sum(shell(u,av,l)*shell(v,av^t,need-l))) for l in (4,6,8))
    if need==12:
        for z in range(4096):
            n+=8*int(np.sum(S[z,np.bitwise_xor.outer(lists[z^a],lists[z^b])^t]))
    return n
expected={5:(3837532,3837340),6:(4943844,4943268),7:(6605112,6603888)}
V=par['V0'];U=V[:,:3]/2;full=V.row_join(U).col_join(U.T.row_join(V[:3,:3]))
for rank,target in expected.items():
    v=full[:rank,:rank];iv=v.inv();bound=[isqrt(int(8*iv[i,i])) for i in range(rank)]
    v2=np.array((2*v).tolist(),dtype=np.int64)
    # Exact integers, vectorized in short batches.
    small={}
    for tail in itertools.product(*(range(-q,q+1) for q in bound[3:])):
        head=list(itertools.product(*(range(-q,q+1) for q in bound[:3])))
        z=np.array([h+tail for h in head],dtype=np.int64)
        norm=np.einsum('ij,ij->i',z@v2,z)
        for row,l in zip(z[(0<norm)&(norm<=16)],norm[(0<norm)&(norm<=16)]):
            assert int(l)%2==0
            mask=sum((int(x)&1)<<i for i,x in enumerate(row))
            key=(mask,int(l)//2);small[key]=small.get(key,0)+1
    out=[]
    for rows in (np.array(B0['old_D4']+B0['new'],dtype=np.int64),raw):
        n=0
        for (mask,l),mult in small.items():
            h=8-l
            p=np.array([mask>>j&1 for j in range(rank)])@rows[:rank]%2
            x,y,z=par['syndrome'](p)
            a,b,t=d['coordA'](x),d['coordA'](y),d['coordB'](z)
            dep=d['dep'](a,b,t)
            assert dep>=h
            if dep==h:
                n+=mult*(count(a,b,t,h) if h else 1)
        out.append(n)
    assert tuple(out)==target,(rank,out,target)
    print('PASS exact norm-8 vectors outside Gamma72',72+rank,'old/new:',out,flush=True)
