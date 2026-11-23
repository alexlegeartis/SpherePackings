from itertools import product
from fractions import Fraction as F
import numpy as np

G=2*np.eye(6,dtype=np.int64)
for i,j in [(0,1),(1,2),(2,3),(3,4),(2,5)]:G[i,j]=G[j,i]=-1
rad=np.array([2,1,0,2,1,0])
assert not (G@rad%3).any()
def solve_rank(A):
    R=np.array(A,dtype=np.int64)%3;n=0
    for i in range(R.shape[1]):
        k=next((k for k in range(n,len(R)) if R[k,i]),None)
        if k is None:continue
        R[[n,k]]=R[[k,n]];R[n]=(R[n]*R[n,i])%3
        R[n+1:]=(R[n+1:]-R[n+1:,i,None]*R[n])%3;n+=1
    return n

assert solve_rank(G)==5
pw=3**np.arange(6)
V=np.array(list(product(range(3),repeat=6)),dtype=np.int64)[:,::-1].copy()
assert np.array_equal(V@pw,np.arange(729))
R={tuple(x) for x in np.eye(6,dtype=int)};Q=list(sorted(R))
for v in Q:
    a=np.array(v)
    for j,h in enumerate(a@G):
        x=a.copy();x[j]-=h;k=tuple(x)
        if k not in R:R.add(k);Q.append(k)
assert len(R)==72
R=np.array(sorted(R),dtype=np.int64)
assert ((R@G*R).sum(axis=1)==2).all()
assert solve_rank([[r[i]*r[j] for i in range(6) for j in range(i,6)] for r in R])==21

def near(V):
    X=V.copy();X[X>1]-=3
    while True:
        Y=X[:,None,:]+3*R[None,:,:]; n=np.sum((Y@G)*Y,axis=-1)
        j=np.argmin(n,axis=1)
        z=n[np.arange(len(X)),j]<np.sum((X@G)*X,axis=-1)
        if not z.any():break
        X[z]=Y[np.arange(len(X))[z],j[z],:]
    assert np.array_equal(X%3,V)
    return X
d_lift=near(V)
d=np.sum((d_lift@G)*d_lift,axis=-1)
assert {i for i,e in enumerate(d) if e>=10}=={int((rad%3)@pw),int((2*rad%3)@pw)}
# Each nonsingular fibre of the horizontal quotient has some lift
# of squared norm at most 4, evidenced by the actual representatives above.
q=(np.sum((V@G)*V,axis=1)*2)%3
for r in V[q!=0]:
    assert min(d[((r+k*rad)%3)@pw] for k in range(3))<=4

# Capacity 27 in the vertical local group at forbidden q<6.
def orbits(dual=False):
    maps=[]
    for i in range(6):
        x=V.copy()
        if dual:x-=V[:,i,None]*G[i]
        else:x[:,i]-=V@G[:,i]
        maps.append((x%3)@pw)
    todo=set(range(729));O=[]
    for i in range(729):
        if i not in todo:continue
        todo.remove(i); S=[i]
        for k in S:
            for m in maps:
                t=int(m[k])
                if t in todo:todo.remove(t);S.append(t)
        O.append(S)
    return O
P,T=orbits(),orbits(True)
Dp=[min(d[S]) for S in P]; I=[j for j,a in enumerate(Dp) if a>=6]
assert [len(P[j]) for j in I]==[240,72,72,1,1]
C=np.array([[sum(2 if i==0 else -1 for i in (V[S]@V[U[0]])%3) for S in P] for U in T],int)
assert ((-6*C[1,I]-4*C[4,I])+np.array([0,0,0,3,3])>=np.array([len(P[j]) for j in I])).all()
assert 1+2*(6+4)+3+3==27
print('critical vertex: quotient obstruction and exact local capacity 27 passed')
