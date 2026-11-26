"""Exact finite-field and integer-HNF replay: every nontrivial character flex of either CHLWW triangular bank is impossible, for all periods. Requires NumPy, SciPy and SymPy. No numerical feasibility thresholds."""
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


from sympy import Matrix
import os, json,sys

import scipy.linalg as la
B=np.array([[21,9,9,12,12],[9,29,13,20,16],[9,13,29,16,20],[12,20,16,32,16],[12,16,20,16,32]])
W=np.array([[1,-1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,1,-1,1,1,3,-1,-1,-1,1,-1],[-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1,-1,-1,-1,3,-1,1,-1,-1,1,-1],[-1,1,-3,-1,-1,1,1,-1,-1,1,-1,-1,1,-1,-1,-1,1,-1,-1,1,-1,-1,-1,-1],[0,2,0,2,2,0,0,0,0,0,0,0,0,0,2,0,0,-2,0,-2,0,2,0,2],[1,-3,1,1,-1,-1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,-1,1]],dtype=int)
H0=480*np.eye(24,dtype=int)-W.T@B@W
sel=list(Matrix(H0).rref()[1]);assert len(sel)==19
a,b=np.triu_indices(19);C=H0[sel][:,sel][a,b]
h=la.null_space(W).T
labs=np.vstack((np.zeros(5,dtype=int),np.eye(5,dtype=int)));P=roots@W.T//8

from math import gcd
from itertools import combinations
from scipy.linalg import qr
p=1009 # prime, used only to prove nonvanishing/no 2-level solution
def pivs(A,p):
 A=A.copy()%p;n,m=A.shape;r=0;pcs=[]
 for j in range(m):
  u=np.flatnonzero(A[r:,j])
  if not len(u):continue
  idx=r+u[0];A[[r,idx]]=A[[idx,r]]
  row=A[r]*pow(int(A[r,j]),-1,p)%p
  for st in range(r+1,n,200):
   ed=min(n,st+200);A[st:ed,j:]=(A[st:ed,j:]-A[st:ed,j:j+1]*row[j:])%p
  r+=1;pcs.append(j)
  if r==n:break
 return pcs
def invp(A):
 n=len(A);R=np.hstack([A.astype(np.int64),np.eye(n,dtype=np.int64)])%p
 for i in range(n):
  u=i+int(np.flatnonzero(R[i:,i])[0]);R[[u,i]]=R[[i,u]]
  R[i]=(R[i]*pow(int(R[i,i]),-1,p))%p
  f=R[:,i].copy();f[i]=0;R=(R-f[:,None]*R[i])%p
 assert np.array_equal(R[:,:n],np.eye(n,dtype=int))
 return R[:,n:]
sg=(1-2*((np.arange(2**17)[:,None]>>np.arange(18))&1)).astype(np.int16) # last sign fixed to +1
rg=np.random.default_rng(26110407)
def good(Z):
 # Z=120(x)+(60)d-2(qB)W in ambient coordinates; 2*60*sqrt(8) times physical midpoint
 if len(Z)<19:return None
 cols=pivs(Z,p); 
 if len(cols)!=18:return None
 Q=Z[:,cols];rr=pivs(Q.T,p);assert len(rr)==18
 cc=(Q@invp(Q[rr]))%p
 test=np.array([k for k in range(len(Q)) if k not in rr],dtype=int)
 # pick strongest constraint first
 test=test[np.argsort(np.count_nonzero(cc[test],axis=1))[::-1]]
 keep=np.arange(2**17)
 for i in test:
  # streaming Gray expansion of dot with all sign choices, without allocating a large matrix product
  if len(keep)>2000:
   A=np.zeros(2**17,dtype=np.int32);A[0]=cc[i].sum()
   for k,c in enumerate(cc[i,:17]):A[1<<k:2<<k]=A[:1<<k]-2*int(c)
   v=A[keep]%p
  else:v=sg[keep]@cc[i]%p
  keep=keep[(v==1)|(v==p-1)]
  if not len(keep):break
 return (cols,rr,keep)

from sympy.matrices.normalforms import hermite_normal_form
normal0=(P==0).all(1);LF=roots[normal0].astype(np.int64);discard=sorted(set(range(24))-set(sel))
target=abs(int(Matrix(W[:,discard].tolist()).det()))*8**7
def index(A): # lattice index for unscaled vectors nominally in L
 H=hermite_normal_form(Matrix(A.T.tolist())); 
 if H.shape!=(19,19):return 0
 assert H.is_upper and target>0
 vols=abs(int(np.prod(np.array(H.diagonal().tolist()[0],dtype=object))));ix,rem=divmod(vols,target);assert rem==0
 return ix
def normprojective(b):return tuple(b*pow(int(b[np.flatnonzero(b)[0]]),-1,p)%p)
for ti,tri in enumerate([[0,4,5],[1,2,3]]):
 E=[];Y=[];certs=[]
 for i,j in combinations(tri,2):
  q=labs[j]-labs[i];xs=roots[(P==q).all(1)].astype(np.int64);z=60*xs-(q@B)@W;ZZ=z[:,sel]
  assert len(pivs((ZZ[:,a]*ZZ[:,b]).T,p))==190
  E.append(xs);Y.append(ZZ)
 potential=set();all_total=0
 for ei,(i,j) in enumerate(combinations(tri,2)):
  q=labs[j]-labs[i];xs=E[ei];z=60*xs-(q@B)@W;sds=[];gd=0;nm=0
  for idx in rg.permutation(len(LF)):
   d=LF[idx];Z=2*z[xs@d==-16]+60*d;re=good(Z)
   if re is None:continue
   cols,rr,kp=re;sds.append(int(idx))
   if len(kp):
    mat=np.vstack((Z[rr][:,sel],60*d[sel]));inv=invp(mat); 
    both=np.vstack((sg[kp],-sg[kp]));bvals=np.hstack((both,np.ones((len(both),1),dtype=int)))@inv.T%p
    for bb in bvals:potential.add(normprojective(bb))
    nm+=2*len(kp)
   if len(sds)>=19:
    gd=index(LF[sds][:,sel])
    if gd==1:break
  print('CAND',ti,ei,'affine generator roots',len(sds),'index',gd,'signatures',nm,'unique cumul',len(potential),flush=True);assert gd==1
 # For each projective candidate: supports across all three edges have enough differences to generate full L?
 bad=[];maxnn=0;maxix=0
 for bcode in sorted(potential):
  b1=np.array(bcode);pool=[];sup=[]
  for xs,zz in zip(E,Y):
   support=xs[(zz@b1%p)!=0][:,sel];sup.append(len(support));pool.append(support[1:]-support[0])
  pool=np.vstack(pool);perm=rg.permutation(len(pool));gd=0
  for nn in [36,65,100,200,400]:
   gd=index(pool[perm[:nn]])
   if gd==1:break
  if gd!=1:bad.append([bcode,gd,sup])
  maxnn=max(nn,maxnn)
 print('EXACT-BLOCH',ti,'tested',len(potential),'nontrivial support groups',len(bad),'nmax',maxnn,'SAMPLES',bad[:10],flush=True)
 assert not bad
