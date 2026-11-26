"""Verify fixed 18D discriminant carrier and every index-two thinning. Requires Numpy and Sympy."""
import numpy as np,itertools
from sympy import Matrix
rows=[]
for i,j in itertools.combinations(range(24),2):
 for a,b in itertools.product((-4,4),repeat=2):
  x=[0]*24;x[i]=a;x[j]=b;rows.append(x)
C=[]
for n in range(4096):
 c=0
 for j in range(12):
  if n>>j&1:c^=0xAE3<<j
 c|=(c.bit_count()%2)<<23;C.append(c)
 if c.bit_count()==8:
  S=[i for i in range(24) if c>>i&1]
  for s in range(256):
   if s.bit_count()%2:continue
   x=[0]*24
   for j,i in enumerate(S):x[i]=2-4*(s>>j&1)
   rows.append(x)
 for j in range(24):
  x=[1-2*(c>>i&1) for i in range(24)];x[j]=-3*x[j];rows.append(x)
roots=np.array(rows,dtype=np.int16);del rows
assert len(roots)==196560

W=np.array([[1,1,1,-1,-1,1,1,1,-1,3,1,1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,-1],[-1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,1,1,1,-1,-1,3,1,1,-1,1],[0,2,2,0,0,0,0,-2,0,2,0,0,0,0,0,0,2,0,-2,0,0,0,-2,-2],[0,2,0,0,2,0,0,0,0,0,0,0,0,-2,0,-2,2,2,-2,0,0,0,0,-2],[2,0,0,0,-2,0,0,0,0,2,-2,0,2,2,0,0,0,0,0,0,0,-2,-2,0],[-1,1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,-1,3,1,1,1,-1]],dtype=np.int32)
claimed=np.array([[0,0,0,0,0,0],[0,81,0,0,0,162],[81,81,0,0,0,0],[0,0,81,0,0,162],[81,0,0,162,0,0],[81,81,162,81,0,162]],dtype=int)
assert roots.shape==(196560,24) and len(C)==len(set(C))==4096 and sum(x.bit_count()==8 for x in C)==759
roots=np.ascontiguousarray(roots,dtype=np.int16)
key={tuple(x):i for i,x in enumerate(roots)}
assert len(key)==196560 and all(tuple(x) in key for x in W)
G=W@W.T//8; M=Matrix(G.tolist());D=int(M.det());assert D==243
J=np.array(M.adjugate().tolist(),dtype=np.int32)
assert np.array_equal(G@J,D*np.eye(6,dtype=int))
P=roots@W.T//8; tw=(16**np.arange(6,dtype=np.int64)); val=(P+4)@tw
U,idx=np.unique(val,return_inverse=True);lookup={int(a):i for i,a in enumerate(U)}
assert all(v in lookup for v in (np.eye(6,dtype=int)+4)@tw) # primitive and all labels realized
Q=np.indices((7,)*6,dtype=np.int16).reshape(6,-1).T-3; prod=Q@J; eng=(prod*Q).sum(1);order=np.argsort(eng)
V,k=np.unique(prod[order]%D,axis=0,return_index=True)
Y=Q[order[k]]; w=eng[order[k]]; assert len(V)==243 and max(w)==2*D
assert max(w)<4*D # every possibly shorter v has |v_i|<=sqrt(2*4)<3 by Cauchy, inside [-3,3]
idxV={tuple(x):i for i,x in enumerate(V)};diff=lambda a,b:idxV[tuple((V[a]-V[b])%D)]
assert list(zip(*[r.tolist() for r in np.unique(w,return_counts=True)]))==[(0,1),(162,72),(324,90),(486,80)]
def echelon(B,v):
 v=int(v)
 while v:
  j=v.bit_length()-1
  if j in B:v^=B[j]
  else:B[j]=v;return True
 return False
rg=np.random.default_rng(91);pr=roots[rg.choice(len(roots),256,replace=False)];char=[];E={}
for z in roots[rg.choice(len(roots),400,replace=False)]:
 p=((pr@z//8)%2).astype('uint8')
 if echelon(E,int.from_bytes(np.packbits(p,bitorder='little').tobytes(),'little')):char.append(z)
 if len(char)==24:break
assert len(char)==24
cm=(roots@np.array(char).T//8%2).astype(np.int64)@(1<<np.arange(24,dtype=np.int64))
ranks=[]
for v in Y:
 m=np.flatnonzero(idx==lookup[int((v+4)@tw)]); assert len(m)
 b={};f=int(cm[m[0]])
 for ii in m[1:]:echelon(b,int(cm[ii])^f)
 assert len(b)==18;ranks.append(len(b))
# All distinct-class distances occur in BOTH halves of every index-two horizontal refinement.
# For zero class the norm-4 section has rank-18 affine residue differences too, so in particular it meets both parities.
C0=idxV[(0,)*6];two3=[i for i in range(D) if w[i]==162];
N=[sum(1<<j for j,b in enumerate(two3) if j!=a and w[diff(v,b)]==162) for a,v in enumerate(two3)]
def has(P,k):
 if k==0:return True
 while P.bit_count()>=k:
  low=P&-P; i=low.bit_length()-1;P^=low
  if has(P&N[i],k-1):return True
 return False
assert not has((1<<len(two3))-1,6)
pts=[idxV[tuple(a)] for a in claimed];assert len(set(pts))==6 and all(w[diff(a,b)]<=162 for a in pts for b in pts)
T=[i for i in range(D) if w[i]==486];plane=None
for a,b in itertools.combinations(T,2):
 W9={idxV[tuple((i*V[a]+j*V[b])%D)] for i,j in itertools.product(range(3),repeat=2)}
 if len(W9)==9 and W9-{C0}<=set(T) and all(np.all(3*V[t]%D==0) for t in W9):plane=W9;break
assert plane is not None # at beta<2 at most 243/9 = 27 points 
# Compare squares exactly: k^2 (1-beta/4)^18 / 243, common denominator.
from fractions import Fraction as F
terms=[1,6**2*F(5,6)**18,27**2*F(2,3)**18,243**2*F(1,2)**18]
assert max(terms)==terms[1]
print("verified 243 nearest classes; affine parity ranks",set(ranks),"baseline crit clique=6, isotropic plane size",len(plane),"exact optimum 6*(5/6)^9/sqrt(243)")
