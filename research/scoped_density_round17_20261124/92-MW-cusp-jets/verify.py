#!/usr/bin/env python3
"""Integer finite arithmetic replay. NumPy is used only for bounded exact integers."""
from pathlib import Path
import json,gzip,numpy as np
p=47
D=json.loads(gzip.decompress(Path(__file__).with_name('input.json.gz').read_bytes()))
G,W,T,U,C,A=[np.array(D[k],dtype='i8') for k in ('G','W','T','U','C','A')]
n=92; I=np.eye(n,dtype='i8')
def mm(X,Y):
    assert X.shape[-1]==Y.shape[0]
    assert X.shape[-1]*int(np.max(np.abs(X),initial=0))*int(np.max(np.abs(Y),initial=0))<2**62
    return X@Y
def eq(X,Y):
    assert np.array_equal(X,Y)
def rank(X,p):
    X=X.copy()%p;r=0
    for j in range(X.shape[1]):
        s=next((i for i in range(r,len(X)) if X[i,j]),None)
        if s is None:continue
        X[[r,s]]=X[[s,r]];X[r]=X[r]*pow(int(X[r,j]),-1,p)%p
        X[r+1:]=(X[r+1:]-X[r+1:,j,None]*X[r])%p;r+=1
    return r
def bareiss(A):
    a=[[int(x) for x in r]for r in A];c=1
    for k in range(len(a)-1):
        b=a[k][k];assert b>0
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                u,v=divmod(b*a[i][j]-a[i][k]*a[k][j],c);assert not v
                a[i][j]=u
        c=b
    return a[-1][-1]
eq(G,G.T); assert bareiss(G)==47**14
assert rank(W,101)==92 and rank(A,p)==2
t=I.copy(); H=np.zeros_like(I); S=np.zeros_like(I)
for i in range(p):
    H-=i*t; S+=t; t=mm(t,T)
eq(t,I);eq(S,np.zeros_like(S))
E=I-T
eq(mm(E,H),p*I);eq(mm(T,A)%p,A)
eq(mm(mm(T,G),T.T),G)
# Source geometrical verification. a+47b encodes a+bi, i^2=-1.
def ae(a,b):return a%p+(b%p)*p
def aplus(a,b):return ae(a%p+b%p,a//p+b//p)
def amin(a,b):return ae(a%p-b%p,a//p-b//p)
def times(a,b):
    x,y=divmod(a,p);z,w=divmod(b,p)
    return ae(y*w-x*z,y*z+x*w)
def apow(a,e):
    r=1
    while e:
        if e&1:r=times(r,a)
        a=times(a,a);e//=2
    return r
inv=[0]+[apow(a,p*p-2) for a in range(1,p*p)]
def trim(P):
    while P and not P[-1]:P.pop()
    return P
def add(P,Q):
    return trim([aplus(P[i] if i<len(P) else 0,Q[i] if i<len(Q) else 0)
                 for i in range(max(len(P),len(Q)))])
def sub(P,Q):return add(P,[amin(0,x)for x in Q])
def mul(P,Q):
    if not P or not Q:return []
    R=[0]*(len(P)+len(Q)-1)
    for i,a in enumerate(P):
        for j,b in enumerate(Q):
            R[i+j]=aplus(R[i+j],times(a,b))
    return trim(R)
def power(P,k):
    R=[1]
    while k:
        if k&1:R=mul(R,P)
        P=mul(P,P);k//=2
    return R
def scale(P,c):return trim([times(x,c) for x in P])
def mod(P,Q):
    R=P[:];c=inv[Q[-1]]
    while len(R)>=len(Q):
        d=len(R)-len(Q);t=times(c,R[-1])
        for i,u in enumerate(Q):R[d+i]=amin(R[d+i],times(t,u))
        trim(R)
    return R
def gcddeg(P,Q):
    while Q:P,Q=Q,mod(P,Q)
    return len(P)-1
def shift(P):
    R=[]
    for c in P[::-1]:R=add(mul(R,[1,1]),[c])
    return R
f=[0]*48;f[1]=46;f[47]=1
i=p; t0=D['cusp']; assert apow(t0,48)==46
k=next(a for a in range(1,p*p) if apow(a,6)==times(45,i))
sections=[];rho=[]
for a,b,c,d in D['unitary']:
    x=[amin(times(a,times(t0,i)),times(b,i)),aplus(times(a,t0),b)]
    y=[amin(times(c,times(t0,i)),times(d,i)),aplus(times(c,t0),d)]
    X=scale(power(trim(x),16),amin(0,apow(inv[k],2)))
    Y=scale(power(trim(y),24),apow(inv[k],3))
    assert sub(mul(Y,Y),power(X,3))==f
    assert len(X)==17 and len(Y)==25 # smooth nonzero specializations
    sections.append((X,Y))
    # Old chart x/y parameter, as in the frozen additive map.
    xr=amin(0,apow(x[1],16)); yr=apow(y[1],24)
    r=times(xr,inv[yr]);rho.append((r%p,r//p))
eq(mm(W,A)%p,np.array(rho,dtype='i8'))
def ht(s,t):
    a,b=s;c,d=t
    x,y=sub(a,c),sub(b,d)
    if not x and not y:return 16
    if not x and not add(b,d):return -16
    h=min([17-len(x)]if x else [1000] ,
          [25-len(y)]if y else [1000])[0]
    return 8-gcddeg(x,y)-h
M=mm(mm(W,G),W.T)
assert all(ht(s,sections[j])==int(M[i,j]) for i,s in enumerate(sections) for j in range(i+1))
M=mm(mm(mm(W,T),G),W.T)
for i,(x,y) in enumerate(sections):
    s=(shift(x),shift(y))
    assert all(ht(s,t)==int(M[i,j]) for j,t in enumerate(sections))
print('Exact Gram, cusp reduction and translation provenance verified')
J=[]
for z in D['first'].values():
    x,y=z['X'],z['Y']
    assert len(trim(x[:]))==19 and len(trim(y[:]))==28
    assert sub(mul(y,y),power(x,3))==f
    a=times(amin(0,x[-1]),inv[y[-1]]);J.append((a%p,a//p))
assert rank(np.array(J),p)==2
def q(v):
    return int(mm(mm(v,G),v.reshape(-1,1))[0])
for depth,norm in [(2,20),(4,24)]:
    v=np.array(D[f'depth{depth}'],dtype='i8')
    assert q(v)==norm
    for i in range(depth):
        v=mm(v,H);assert not np.any(v%p);v=v//p
print('Depth two and four norm witnesses verified; independent first jets verified')
# Just two even-determinant comparisons suffice for every between-level step.
for min_upper,exponent in [(18,2),(20,3),(20,4),(24,5),(24,6),(24,7),(24,8)]:
    assert min_upper**46<16**46*47**exponent
print('All consecutive intermediate sublattices through depth four are below MW92')
# The two-jet quotient as pairs (u,v), with old cusp coordinates u from A.
F=np.array(D['F'],dtype='i8')
i,j=next((i,j)for i in range(n)for j in range(i) if (A[i,0]*A[j,1]-A[i,1]*A[j,0])%p)
a,b,c,d=(int(x) for x in A[[i,j]].ravel())
R=np.array([[d,-b],[-c,a]],dtype='i8')*pow((a*d-b*c)%p,-1,p)@I[[i,j]]%p
eq(mm(R,A)%p,np.eye(2,dtype='i8'))
B=mm(I-mm(A,R),H)
assert not np.any(B%p)
B=mm(B//p,A)%p
J=np.column_stack((A,B))
assert rank(J,p)==4
eq(mm(E,B)%p,A);eq(mm(mm(E,E),J)%p,np.zeros_like(J))
def quotient(X):
    Y=J.copy();M=I.copy()
    for k in range(4):
        t=next(i for i in range(k,n) if Y[i,k])
        Y[[k,t]]=Y[[t,k]];M[[k,t]]=M[[t,k]]
        c=pow(int(Y[k,k]),-1,p);Y[k]=c*Y[k]%p;M[k]=c*M[k]%p
        for t in range(n):
            if t!=k:
                c=Y[t,k];Y[t]=(Y[t]-c*Y[k])%p;M[t]=(M[t]-c*M[k])%p
    Q=mm(mm(M[:4],X),J)%p
    eq(mm(X,J)%p,mm(J,Q)%p)
    return Q
for X in [U,C,F]:eq(mm(mm(X,G),X.T),G)
S,C0,F0=[quotient(X) for X in [U,C,F]]
assert all(not np.any(X[2:,:2]) for X in [S,C0,F0])
QQ=[[int(x)for x in X.ravel()]for X in [S,C0,F0]]
def case(x):
    a,b,c,d=map(int,x)
    if a: return (1,b*pow(a,-1,p)%p,(a*d-b*c)*pow(a*a%p,-1,p)%p)
    assert b
    return (0,1,(a*d-b*c)*pow(b*b%p,-1,p)%p)
Z=np.zeros((p*p,p),dtype=bool)
U16=np.array(D['q16rows'],dtype='i8')
assert len(U16)==17 and all(q(v)==16 for v in U16)
X=mm(U16,J)%p
t=X.copy()
for j in range(92):
    r=t.copy()
    for k in range(3):
        a,b,c,d=r.T;Z[a*p+b,(a*d-b*c)%p]=True
        r=mm(r,C0)%p
    eq(r,t)
    t=mm(t,S)%p
eq(t,X)
assert (Z.any(axis=1)[1:]).all()
def block_choices():
    for a,b in [(1,j)for j in range(p)]+[(0,1)]:
        seen=set()
        for t in range(1,p):
            u,v=t*a%p,t*b%p
            seen.update(int(c)*pow(t*t%p,-1,p)%p for c in np.flatnonzero(Z[u*p+v]))
        for s in seen:yield (a,b,s)
bad=set(block_choices())
X=np.array(D['q18rows'],dtype='i8')
assert all(q(v)==18 for v in X)
P=mm(X,J)%p
assert not np.any(P[:,:2])
Y=P.copy();u18=set()
for j in range(92):
    v=Y.copy()
    for k in range(3):
        u18.update(case([a,b,0,0])[:2] for a,b in v[:,2:])
        v=mm(v,C0)%p
    Y=mm(Y,S)%p
assert len(u18)==6
open_cases={(a,b,c) for a,b in [(1,j) for j in range(p)]+[(0,1)] if (a,b) not in u18
             for c in range(p) if (a,b,c) not in bad}
assert len(bad)==1386 and len(open_cases)==864
# In each residual conjugacy class produce one representative; some extra
# isometries may also send a residual case into the already covered set.
rem=set(open_cases);orbit_sizes=[];orbit_reps=[];closed=0
for a,b,c in sorted(open_cases):
    if (a,b,c) not in rem:continue
    row=np.array([a,b,0,c] if a else [a,b,-c,0],dtype='i8')
    orbit={(a,b,c)};front=[(a,b,c)]
    while front:
        a0,b0,c0=front.pop()
        row=np.array([a0,b0,0,c0] if a0 else [a0,b0,-c0,0],dtype='i8')
        assert case(row)==(a0,b0,c0)
        for X in [S,C0,F0]:
            rv=mm(row,X)%p
            w=mm(np.array([0,0,a0,b0],dtype='i8'),X)%p
            # Check that this map actually normalizes the square-zero
            # translation action up to a scalar for the present plane.
            assert not np.any(w[:2]) and (w[2]*rv[1]-w[3]*rv[0])%p==0
            k=case(rv)
            if k not in orbit:orbit.add(k);front.append(k)
    rem-=orbit
    if orbit<=open_cases:orbit_sizes.append(len(orbit));orbit_reps.append(((a,b,c),len(orbit)))
    else:closed+=len(orbit&open_cases)
print('Remaining translated-line cases: before/after finite isometries',
      len(open_cases),len(open_cases)-closed,'in orbits',sorted(orbit_sizes))
print('Representatives (a,b,c), orbit size:',orbit_reps)
n0=103823
Q=np.indices((p,p,p),dtype='i2').reshape(3,-1).T
Vs=np.array([(a*a%p,a*b%p,b*b%p)for a in range(p)for b in range(p)if a or b],dtype='i8')
FZ=[[int(d)for d in np.flatnonzero(r)]for r in Z[1:]]
def graph_minimum():
    best=2208
    for i in range(0,n0,1200):
        M=mm(Q[i:i+1200].astype('i8'),Vs.T)%p
        count=np.zeros(len(M),dtype='i8')
        for j,L in enumerate(FZ):
            for z in L:count+=(M[:,j]==z)
        best=min(best,int(np.min(count)))
    return best
assert len(Q)==n0==p**3
assert graph_minimum()==68
print('All 103823 graph quadratics hit >=68 literal norm-16 residues; q16 hits every nonzero u')
