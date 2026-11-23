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
    maps.append(((-V)%q)@pw)
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
    return np.sum((a@G)*a,axis=-1),np.sum((b@G)*b,axis=-1),np.sum((a@G)*b,axis=-1)
anchors=V[[s[0] for s in P]].reshape(len(P),2,n).astype(np.int64)
X=anchors.copy();X[X>q/2]-=q
lines=[set() for _ in P]
for c in [16,18,20,24,28,32,36,40,48,56,64,80,96]:
 for t in [0,2,4,6,8,10,12,14,16]:
    x=X.copy()
    while True:
        a,b,d=parts(x); score=16*a+c*b+t*d
        y=x[:,None,:,:]+q*R[None,:,:,:]
        u,v,w=parts(y); vals=16*u+c*v+t*w
        ix=np.argmin(vals,axis=1)
        z=vals[np.arange(len(P)),ix]<score
        if not z.any():break
        x[z]=y[np.arange(len(P))[z],ix[z],:,:]
    assert np.array_equal(x%q,anchors)
    a,b,d=parts(x)
    for s,z in zip(lines,zip(map(int,a),map(int,b),map(int,d))):s.add(z)
assert lines[0]=={(0,0,0)}

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
if q==3:
    from sympy import Matrix
    from math import isqrt
    from functools import lru_cache
    B=Matrix(G.tolist()).inv()
    @lru_cache(maxsize=None)
    def shell(r, hi):
        lim=[isqrt(int(hi*B[i,i])) for i in range(6)]
        W=np.array(list(it.product(*[[j for j in range(-l,l+1) if j%q==r[i]] for i,l in enumerate(lim)])),dtype=np.int64)
        if not len(W):return W.reshape(0,6),np.array([],dtype=np.int64)
        d=np.sum((W@G)*W,axis=-1)
        return W[d<=hi],d[d<=hi]
    focus=[(F(9,4),F(9,4),F(0)),(F(3,2),F(3,2),F(3,2)),(F(9,5),F(18,5),F(0)),
        (F(36,25),F(81,50),F(36,25)),(F(18,11),F(18,11),F(18,11)),(F(3,2),F(11,7),F(10,7))]
    for a,c,t in focus:
        import math
        z=math.lcm(a.denominator,c.denominator,t.denominator)
        boundU=int(F(18)/(a-t*t/(4*c)));boundV=int(F(18)/(c-t*t/(4*a)))
        count=0
        for j,S in enumerate(lines):
            if j==0 or any(a*d+c*e+t*f<18 for d,e,f in S):continue
            x,y=anchors[j]
            X,u=shell(tuple(x),boundU);Y,v=shell(tuple(y),boundV)
            if not(len(X) and len(Y)):continue
            C=X@G@Y.T; K=int(a*z)*u[:,None]+int(c*z)*v[None,:]+int(t*z)*C
            i,k=np.unravel_index(K.argmin(),K.shape)
            if K[i,k]<18*z:
                S.add((int(u[i]),int(v[k]),int(C[i,k])));count+=1
        print('exact extra removals',(a,c,t),count,flush=True)

# X=(a-t,c-a,t), metric [[a,t/2],[t/2,c]];
# X>=0, a>=1. Remove inequalities universally satisfied on this cone.
lim=2*q*q
L=[sorted({(a+b,b,a+b+d) for a,b,d in S if min(a+b,a+b+d)<lim}) for S in lines]
L[0]=[]
W=sorted({z for S in L[1:] for z in S})
assert (0,0,0) not in W
lookup={z:i for i,z in enumerate(W)}
ids=[np.array([lookup[z] for z in S],dtype=int) for S in L]
A=[(*r,lim) for r in W]+[(1,0,0,0),(0,1,0,0),(0,0,1,0),(1,0,1,1)]
def det(u,v,w):
    a,b,c=u; d,e,f=v;g,h,i=w
    return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
def verts():
    uniq=set()
    for B in it.combinations(A,3):
        C=list(zip(*B)); E=C[:3]; h=det(*E)
        if not h:continue
        x,y,z=[det(*(E[:j]+C[3:]+E[j+1:])) for j in range(3)]
        if h<0:h,x,y,z=-h,-x,-y,-z
        if min(x,y,z)<0 or x+z<h:continue
        import math
        gcd=math.gcd(math.gcd(x,y),math.gcd(z,h))
        uniq.add((x//gcd,y//gcd,z//gcd,h//gcd))
    return uniq
P0=verts(); print('gates',len(W),'vertices',len(P0),flush=True)
Wn=np.array(W,dtype=np.int64)
masks={}
for x,y,z,h in P0:
    a=x+z;c=x+y+z
    cdet=F(4*a*c-z*z,4*h*h)
    mask=Wn@np.array([x,y,z],dtype=np.int64)>=lim*h
    S=tuple(k for k,sz in enumerate(ids) if k and (mask[sz]).all())
    if S not in masks or (cdet,(F(a,h),F(c,h),F(z,h)))<masks[S]:masks[S]=(cdet,(F(a,h),F(c,h),F(z,h)))
print('masks',len(masks),flush=True)
add=[]
out=[]
for S,(d,X) in masks.items():
    # Cheap cardinality first, exact.
    target=F(64,9)*d**3
    if 1+sum(map(int,sz[list(S)]))<=target:continue
    c=cap(S); k=c.numerator//c.denominator
    out.append((F(k)/target,X,c,len(S)))
    rank=0
    while q**rank<=target:rank+=1
    if q**rank<=k:add.append((rank,X))
    # This control is the explicit four-point binary E6 x E6 code.
    if q==2 and X==(F(1),F(1),F(0)):assert k>=4
    # A code in this graph must obey the critical orthogonal depth gate.
    if q==3 and X==(F(9,5),F(18,5),F(0)):
        from verify_additive_exception import G as G0,d as d0,pw as pw0
        assert (G0==G).all()
        assert all(int(d0[V[P[j][0],:n]@pw0])+2*int(d0[V[P[j][0],n:]@pw0])>=10 for j in S)
print('worst',sorted(out,reverse=True)[:10],flush=True)
print('open',[(str(r),tuple(map(str,t)),float(v),s) for r,t,v,s in sorted(out,reverse=True) if r>1],flush=True)

print('additive candidates',add,flush=True)
if q==2:assert not add and all(r<=1 for r,*_ in out)
if q==3:assert set(add)=={(7,(F(9,5),F(18,5),F(0)))}
