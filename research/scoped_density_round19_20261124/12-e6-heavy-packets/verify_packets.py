"""Exact orthogonal-stabilizer three-point certificate for 12D heavy packets."""
from itertools import product, permutations
from fractions import Fraction as F
from pathlib import Path
from collections import defaultdict
import json
import numpy as np

data=json.loads(Path(__file__).with_name('packet-dual.json').read_text())
N=data['N']; assert N==109
# Coordinates (0,w) in E6/3E6 represent the five-dimensional quotient
# by the radical (2,1,0,2,1,0), with quadratic value Q(x)/2.
coords=list(product(range(3),repeat=5))
V=np.array(coords,dtype=np.int64)
I={x:i for i,x in enumerate(coords)}
E=((0,1),(1,2),(2,3),(3,4),(2,5))
def relation(B):
    X=np.concatenate((np.zeros(B.shape[:-1]+(1,),dtype=np.int64),B),axis=-1)
    q=((X**2).sum(-1)-sum(X[...,i]*X[...,j] for i,j in E))%3
    return np.where(q!=0,q,np.where((B==0).all(-1),0,3))
T=relation(V)
D=relation((V[:,None,:]-V[None,:,:])%3)
P=np.empty((243,243,4),dtype=np.int64)
P[:,:,0]=T[:,None]; P[:,:,1]=T[None,:]; P[:,:,2]=D
P[:,:,3]=((V[:,None,:]+V[None,:,:])%3==0).all(-1)
orbitals, ids=np.unique(P.reshape((-1,4)),return_inverse=True,axis=0)
ids=ids.reshape((243,243))
B=np.array([ids==j for j in range(len(orbitals))],dtype=np.int64)
assert len(B)==40
mult=B.sum((1,2))
counts=np.bincount(T)
assert tuple(counts)==(1,72,90,80)
def perm_ids(u,v):
    X=V[[0,u,v]]
    return tuple(int(ids[I[tuple((X[p[1]]-X[p[0]])%3)],
                         I[tuple((X[p[2]]-X[p[0]])%3)]])
                 for p in permutations(range(3)))
pm=[]
for m in B:
    vec=[perm_ids(u,v) for u,v in np.argwhere(m)]
    assert len(set(vec))==1
    pm.append(vec[0])

def valid(o,p):
    return o==p==0 or (o>0 and p>0 and o+p>=4)
n=40
Q={}
for i,j in product(range(n),repeat=2):
    if all(valid(a,b)for a,b in zip(orbitals[i,:3],orbitals[j,:3])):
        Q[i,j]=min(min((u,v),(v,u)) for u,v in zip(pm[i],pm[j]))
keys=sorted(set(Q.values()))
assert len(keys)==56
I={p:i for i,p in enumerate(keys)}
X=np.zeros((n,n,len(I)),dtype=np.int64)
a=np.zeros((16,len(I)),dtype=np.int64)
Y=np.zeros_like(a)
for (i,j),k in Q.items():
    p=I[k]
    X[i,j,p]=1
    if orbitals[i,0]==orbitals[j,0]==0:
        a[4*orbitals[i,1]+orbitals[j,1],p]=1
    Y[4*orbitals[i,0]+orbitals[j,0],p]+=mult[i]*mult[j]
G=np.r_[[N*N*a[0]],Y-N*a*np.kron(counts,counts)[:,None]]
comp=np.empty_like(X)
for i,j in product(range(n),repeat=2):
    comp[i,j]=a[4*orbitals[i,2]+orbitals[j,2]]-X[i,j]
l=list(map(F,data['l']))
assert len(l)==len(G) and l[0]*N*N==1
cert=[sum(c*int(q) for c,q in zip(l,s))for s in G.T]
U={int(i):np.array(u,dtype=np.int64) for i,u in data['U'].items()}
# All temporary int64 products have fewer than 243^2 terms and factors <=10000.
assert all(len(u)==243 and int(abs(u).max())<=10000 for u in U.values())
H={i:np.einsum('ri,krs,sj->kij',u,B,u).astype(object)for i,u in U.items()}
Z=X.reshape((n*n,-1)).astype(object)
Zc=comp.reshape((n*n,-1)).astype(object)
for i,j,c,v,w in data['rows']:
    w=F(w)
    assert w>0 and c in (0,1)
    L,J=H[i],H[j]
    A=np.asarray(v,dtype=object).reshape((L.shape[1],J.shape[1]))
    # Exact Rayleigh evaluation on (U_i tensor U_j) vec(A).
    q=((L@A).reshape((n,-1))@np.array([(A@m).ravel()for m in J],dtype=object).T).ravel()
    coef=q@(Zc if c else Z)
    cert=[x+w*int(y) for x,y in zip(cert,coef)]
assert max(cert)<0
print('PASS: 56 exact coefficients strictly negative,',len(data['rows']),'positive Rayleigh cut multipliers; no 109-packet code.')

# A literal quotient two-point barrier. Each entry is per displacement.
# Twice the real character sums of each of the four distance shells.
Wreal=set()
for h in V:
    v=(V@h)%3
    Wreal.add(tuple(3*int(np.sum((T==k)&(v==0)))-int(counts[k])for k in range(4)))
assert Wreal=={(2,-18,18,-2),(2,0,-18,16),(2,18,0,-20),(2,144,180,160)}
p={(0,0):F(1),(1,3):F(1,160),(2,2):F(1,150),
   (3,1):F(1,160),(3,3):F(1,800)}
assert all(valid(i,j) for i,j in p)
assert all(sum(h[i]*j[k]*w for (i,k),w in p.items())>=0 for h,j in product(Wreal,repeat=2))
assert sum(int(counts[i])*int(counts[j])*w for (i,j),w in p.items())==135
print('PASS: exact quotient two-point pseudo-size 135')
