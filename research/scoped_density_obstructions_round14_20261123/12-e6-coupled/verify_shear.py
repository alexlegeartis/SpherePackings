import itertools as it, sys
import numpy as np
from scipy.optimize import linprog
from fractions import Fraction as F

q=int(sys.argv[1])
assert q in (2,3)
n=6
G=np.eye(n,dtype=np.int64)*2
for i,j in [(0,1),(1,2),(2,3),(3,4),(2,5)]:G[i,j]=G[j,i]=-1
N=q**(2*n)
pw=np.array([q**i for i in range(2*n)])
V=((np.arange(N)[:,None]//pw)%q).astype(np.int8)

def orbits(dual):
    maps=[]
    for i in range(n):
        a=V.astype(np.int64)
        if dual:
            a[:,:n]-=a[:,i,None]*G[i]
            a[:,n:]-=a[:,n+i,None]*G[i]
        else:
            a[:,i]-=a[:,:n]@G[:,i]
            a[:,n+i]-=a[:,n:]@G[:,i]
        maps.append((a%q)@pw)
    maps.extend((((V[:,n:]%q)@pw[:n])+((V[:,:n]%q)@pw[n:]),((-V)%q)@pw))
    done=np.zeros(N,bool); out=[]
    for g in range(N):
        if done[g]:continue
        S=[g];done[g]=True
        for s in S:
            for m in maps:
                v=int(m[s])
                if not done[v]:done[v]=True;S.append(v)
        out.append(np.asarray(S,dtype=np.int64))
    return out

P=orbits(False); T=orbits(True); print('orbits',q,len(P),len(T),flush=True)
assert len(P)==len(T)
sz=np.asarray([len(S) for S in P],int)
orb=np.empty(N,dtype=np.int64)
for i,S in enumerate(P):orb[S]=i
lut=[2,-2] if q==2 else [2,-1,-1]
M=[]
for Tt in T:
    a=(V@V[Tt[0]])%q
    M.append(sum(lut[r]*np.bincount(orb[a==r],minlength=len(P)) for r in range(q)))
M=np.asarray(M,np.int64)
assert (M[0]==2*sz).all() and (M[:,0]==2).all()

# Build valid short-lift witnesses; missing lifts only relax the gate graph.
R=[np.zeros(n,dtype=np.int64)]
for i in range(n): R.append(np.eye(n,dtype=np.int64)[i])
R0={tuple(x) for x in R if x.any()}; Todo=list(R0)
for v in Todo:
    z=np.asarray(v)
    for j,h in enumerate(z@G):
        y=z.copy();y[j]-=h;k=tuple(y)
        if k not in R0:R0.add(k);Todo.append(k)
assert len(R0)==72
R=np.asarray([(a,b) for a in sorted(R0) for b in [np.zeros(n,dtype=np.int64)]]+
             [(a,b) for b in sorted(R0) for a in [np.zeros(n,dtype=np.int64)]]+
             [(a,-np.asarray(a)) for a in sorted(R0)],dtype=np.int64)
def parts(x):
    a=x[...,0,:];b=x[...,1,:]
    return np.sum((a@G)*a+(b@G)*b,axis=-1),np.sum((a@G)*b,axis=-1)
anchors=V[[s[0] for s in P]].reshape(len(P),2,n).astype(np.int64)
X=anchors.copy();X[X>q/2]-=q
lines=[set() for _ in P]
for p in range(17):
    x=X.copy()
    while True:
        a,b=parts(x); score=16*a+p*b
        y=x[:,None,:,:]+q*R[None,:,:,:]
        u,v=parts(y); vals=16*u+p*v
        ix=np.argmin(vals,axis=1)
        z=vals[np.arange(len(P)),ix]<score
        if not z.any():break
        x[z]=y[np.arange(len(P))[z],ix[z],:,:]
    assert np.array_equal(x%q,anchors)
    a,b=parts(x)
    for s,z in zip(lines,zip(map(int,a),map(int,b))):s.add(z)
assert lines[0]=={(0,0)}

# Distinct allowed masks are tested at every endpoint and at midpoints.
# At fixed mask, covolume decreases in t; use its greatest t, even if
# the endpoint admits extra orbits.
bound=2*q*q
points=sorted({F(0),F(1)}|{F(bound-a,b) for L in lines for a,b in L if b and 0<F(bound-a,b)<1})
masks={}
for t,edge in [(x,x) for x in points]+[((a+b)/2,b) for a,b in zip(points,points[1:])]:
    S=tuple(j for j,L in enumerate(lines) if j and all(F(a)+t*b>=bound for a,b in L))
    masks[S]=max(edge,masks.get(S,F(0)))
def cap(S):
    if not S:return F(1)
    S=list(S)
    A=np.r_[-np.unique(M[1:,S],axis=0),np.diag(sz[S])]
    b=np.r_[np.full(len(A)-len(S),2),sz[S]]
    z=linprog(-sz[S], A_ub=A,b_ub=b,method='highs')
    assert z.success
    Y=[(i,F(max(0,-y)).limit_denominator(10**5)) for i,y in enumerate(z.ineqlin.marginals) if -y>1e-10]
    residual=sum(max(F(0),F(int(sz[k]))-sum(y*int(A[i,j]) for i,y in Y)) for j,k in enumerate(S))
    return 1+sum(int(b[i])*y for i,y in Y)+residual
def strict_target(t):
    return F(64,9)*(1-t*t/4)**3
res=[]
for S,t in masks.items():
    val=cap(S);c=val.numerator//val.denominator
    res.append((F(c)/strict_target(t),t,val,len(S)))
    if t==1 and q==3:print("terminal witness",[(V[P[k][0]].tolist(),sorted(lines[k])) for k in S],flush=True)
print('points',len(points),'masks',len(masks),'worst',sorted(res,reverse=True)[:12],flush=True)
print('open',[(str(r),str(t),str(v),s) for r,t,v,s in res if r>1],flush=True)
assert all(r<=1 for r,*_ in res)
if q==3:assert max(r for r,*_ in res)==1
