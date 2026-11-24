from itertools import product
from pathlib import Path
import numpy as np
V=np.array(list(product(range(3),repeat=5)),dtype=np.int64)
G=2*np.eye(5,dtype=np.int64)
for i,j in [(0,1),(1,2),(2,3),(1,4)]:G[i,j]=G[j,i]=-1
def typ(D):
    Q=((D@G)*D).sum(-1)//2%3
    return np.where(Q,Q,np.where((D==0).all(-1),0,3))
T=typ((V[:,None]-V[None])%3)
for name,M,p in [('aperiodic-57.txt',57,1),('periodic-63.txt',63,9)]:
    A=np.loadtxt(Path(__file__).with_name(name),dtype=np.int64)
    assert A.shape==(M,2) and (A>=0).all() and (A<243).all()
    x,y=A.T;a=T[x[:,None],x[None,:]];b=T[y[:,None],y[None,:]]
    off=~np.eye(M,dtype=bool)
    assert (a[off]>0).all() and (b[off]>0).all() and ((a+b)[off]>=4).all()
    C=V[A].reshape((-1,10));S=set(map(tuple,C))
    shifts=[v for v in (C-C[0])%3 if all(tuple(c)in S for c in (C+v)%3)]
    assert len(shifts)==p
    K=np.ones((243,243),bool)
    for s,t in A:
        a=T[s,:,None];b=T[t,None,:]
        K&=(a>0)&(b>0)&(a+b>=4)
    assert not K.any()
    print('PASS',M,'labels,',len(shifts),'translations in stabilizer, no addable quotient label,',M*9,'sites')
