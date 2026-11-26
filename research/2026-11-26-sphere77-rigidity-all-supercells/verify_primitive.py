"""Exact verification of 19D CHLWW internal-bank equalities, ranks, and positive stresses. Requires NumPy, SciPy, SymPy; O(2 s) single thread. Floating-point arithmetic only guesses a partition and a rational solution. The proof tests are exact and do not require trusting conditioning."""
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
def rankmod(A,p=32749):
 A=A.copy()%p;n,m=A.shape;r=0
 for j in range(m):
  u=np.flatnonzero(A[r:,j])
  if not len(u):continue
  idx=r+u[0]; A[[r,idx]]=A[[idx,r]]; row=A[r]*pow(int(A[r,j]),-1,p)%p
  for st in range(r+1,n,400): 
   end=min(n,st+400);A[st:end,j:]=(A[st:end,j:]-A[st:end,j:j+1]*row[j:])%p
  r+=1
  if r==n:break
 return r
for tri in [[0,4,5],[1,2,3]]:
 fs=[];zs=[];qs=[];obs=[]
 for i,j in itertools.combinations(range(3),2):
  q=labs[tri[j]]-labs[tri[i]];xs=roots[(P==q).all(1)]
  assert len(xs)==2011
  v=xs@h.T/np.sqrt(8);z=60*xs.astype(np.int64)-(q@B)@W;z=z[:,sel]
  mat=np.zeros((len(xs),228));mat[:,:190]=v[:,a]*v[:,b]*(1+(a!=b))
  raw=np.zeros_like(mat,dtype=np.int64);raw[:,:190]=19*z[:,a]*z[:,b]
  for k,e in [(i,-1),(j,1)]:
   if k:mat[:,190+19*(k-1):190+19*k]=e*v;raw[:,190+19*(k-1):190+19*k]=e*z
  fs.append(mat);zs.append(raw)
  obs.append(xs)
 F=np.concatenate(fs);R=np.concatenate(zs)
 target=np.zeros(228);target[:190]=(a==b)*len(R)*(52/15)/19
 w=F@la.solve(F.T@F,target,assume_a='pos')
 classes,groups=np.unique(np.rint(w*10**7).astype(int),return_inverse=True)
 # certify column rank of the UNcentered contact matrix by an exact finite field calculation
 rk=rankmod(R[:,:].T);assert rk==228
 R[:,:190]-=208*C
 A=np.zeros((len(classes),228),dtype=np.int64);np.add.at(A,groups,R);A=A.T
 Q,R0,piv=la.qr(A.astype(float),pivoting=True,mode='economic')
 dg=np.abs(np.diag(R0));r=np.sum(dg>dg.max()*1e-10);cols=piv[:r]
 _,_,piv2=la.qr(A[:,cols].T.astype(float),pivoting=True,mode='economic')
 ix=piv2[:r];E=Matrix(A[ix][:,cols].tolist()); RR=Matrix(A[ix].tolist());p0=Matrix(classes.tolist());v = E.inv()*(-RR*p0)
 for t,k in enumerate(cols):p0[int(k)]+=v[t]
 assert min(p0)>0 and Matrix(A.tolist())*p0==Matrix.zeros(228,1)
 print('EXACT',tri,'rank',rk,'group rank',r,'classes',len(classes),'min weight',float(min(p0))/1e7,'max digits',max(len(str(t.q)) for t in p0),flush=True)
