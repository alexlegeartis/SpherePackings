"""Exact certificate, 20D and 21D Leech section realizations. Requires NumPy and SymPy."""
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


from fractions import Fraction as F
Cb=np.array([[-2,2,2,0,2,-2,0,2,0,0,0,0,2,0,0,0,0,0,0,0,0,0,2,0],[0,-2,0,0,0,0,2,0,0,0,-2,2,0,0,0,-2,0,0,-2,0,0,0,-2,-2],[-1,1,1,-1,-1,-1,-3,-1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,1],[0,0,2,-2,0,-2,0,0,-2,0,0,-2,0,0,0,0,0,-2,-2,0,2,0,0,0],[-1,-1,-1,3,-1,1,-1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,-1,-1,1,1,1,-1],[-1,-1,-1,-1,1,-1,3,-1,1,-1,-1,1,-1,1,1,-1,1,1,1,1,-1,1,1,1],[-1,-1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,-3,-1,-1,-1,-1,1,-1,1,-1,1,-1],[0,-2,2,-2,0,0,0,0,-2,0,0,0,0,0,-2,-2,0,0,2,0,0,0,0,-2],[1,3,1,-1,1,-1,-1,1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1,1,-1,1,1,1],[0,-2,0,2,-2,-2,-2,0,0,0,0,0,-2,0,0,-2,0,0,0,0,0,0,0,2],[1,1,-1,1,-1,-1,1,1,1,1,-1,-3,-1,1,1,1,1,1,1,-1,1,-1,1,-1],[1,1,-1,-1,1,1,3,1,-1,1,-1,-1,1,1,-1,1,1,1,1,1,1,1,-1,1],[0,0,2,0,-2,-2,0,0,0,0,0,0,0,0,0,-2,0,0,-2,0,2,0,2,2],[0,0,-2,0,2,0,0,0,-2,0,2,0,2,0,0,0,0,2,0,0,2,0,2,0],[1,-1,1,-1,-1,1,1,1,1,1,1,3,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1],[0,-2,0,0,-2,0,-2,0,2,0,0,-2,0,-2,2,0,0,0,0,0,0,0,0,-2],[1,1,-1,-1,1,1,1,3,-1,1,-1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,-1],[-1,1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,-3,1,-1,1,-1,1,1,-1,-1,-1,1],[-2,0,2,0,2,0,0,2,0,0,0,0,0,0,2,2,-2,0,2,0,0,0,0,0],[-1,-1,1,-1,-3,1,1,1,1,-1,1,1,1,-1,-1,-1,1,-1,-1,1,-1,-1,1,-1],[-1,-1,1,-1,1,-1,-1,1,1,3,1,1,-1,-1,1,1,1,1,-1,-1,-1,1,1,-1],[0,0,0,0,2,0,0,0,0,-2,2,0,0,0,0,0,0,2,2,2,0,0,2,-2],[0,0,0,0,-2,0,0,2,2,0,0,2,0,-2,0,-2,-2,0,0,2,0,0,0,0],[-1,1,-1,1,-1,-1,-3,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,1]],dtype=np.int16);known={row.tobytes() for row in roots};assert all(row.tobytes() in known for row in Cb)
assert not np.any((roots@Cb.T)%8)
P=(roots@Cb.T//8%2).astype(np.int64)@(1<<np.arange(24,dtype=np.int64));
def modrank(vs):
 b={}
 for v in vs:
  v=int(v)
  while v:
   j=v.bit_length()
   if j not in b:b[j]=v;break
   v^=b[j]
 return len(b)
assert modrank(P)==24
def maxclq(A):
 best=0;calls=0
 def go(P,k):
  nonlocal best,calls
  calls+=1;best=max(k,best)
  vs=[];bo=[];i=0;l=P
  while l:
   q=l;i+=1
   while q:
    bit=q&-q;j=bit.bit_length()-1;vs.append(j);bo.append(i);q&=~(A[j]|bit);l^=bit
  for j,b in zip(vs[::-1],bo[::-1]):
   if k+b<=best:return
   go(P&A[j],k+1);P^=1<<j
 go((1<<len(A))-1,0);return best,calls

from itertools import combinations
for rank,ids in [(3,[23,168384,140814]),(4,[23,172938,146593,158539])]:
 d=24-rank;W=roots[ids].astype(np.int64);G=W@W.T//8
 H=Matrix(G.tolist());D=int(H.det());J=np.array(H.adjugate().tolist(),dtype=np.int64)
 assert D==({3:36,4:125}[rank])
 P8=8*D*np.eye(24,dtype=np.int64)-W.T@J@W
 cols=list(Matrix(P8.tolist()).rref()[1]);assert len(cols)==d
 raw=roots.astype(np.int64)@W.T;assert not (raw%8).any(); profiles=raw//8 
 sites=np.vstack([np.zeros(rank,dtype=np.int64),np.eye(rank,dtype=np.int64)])
 assert all(np.any(np.all(profiles==q,axis=1)) for q in sites[1:])
 beta=max((a-b)@J@(a-b) for a,b in combinations(sites,2));assert beta==({3:15,4:50}[rank])
 edges=[];done=False
 for i,j in combinations(range(len(sites)),2):
  q=sites[i]-sites[j]
  if q@J@q!=beta:continue
  rr=roots[np.all(profiles==q,axis=1)].astype(np.int64);R=D*rr-(q@J@W)[None,:]
  E=R.T@R;assert np.array_equal(E*(d*8*D),E.trace()*P8) and not np.any(R.sum(0))
  assert len(R)==({3:9072,4:4750}[rank]);edges.append((i,j))
  if not done:
   ix=np.triu_indices(d);m=d*(d+3)//2
   u=R[np.random.default_rng(1927).choice(len(rr),2*m,replace=False)][:,cols]
   A=np.hstack((u[:,ix[0]]*u[:,ix[1]],u))%65521
   for c in range(m):
    piv=np.flatnonzero(A[c:,c]);assert len(piv);h=c+int(piv[0]);A[[c,h]]=A[[h,c]]
    A[c]=A[c]*pow(int(A[c,c]),-1,65521)%65521
    A[c+1:]=(A[c+1:]-A[c+1:,c,None]*A[c])%65521
   done=True
 assert (set(edges)=={(0,1),(0,2),(1,3),(2,3)} if rank==3 else len(edges)==10)
 print(d,"D =",D,", beta =",beta,"/",D,", edges =",edges,", shell size =",len(R),", rank =",m)

 Q=np.indices((9,)*rank).reshape(rank,-1).T-4;V=Q@J;energy=(Q*V).sum(1);o=np.argsort(energy,kind='stable');U,ix=np.unique(V[o]%D,axis=0,return_index=True);w=energy[o[ix]];qmin=Q[o[ix]];assert len(U)==D and max(w)<4*D
 # Outside this box, Cauchy gives normal squared length >=25/4 (all normal roots have squared length 4).
 key=lambda v:(v+4)@9**np.arange(rank)
 vals=key(profiles);mp={tuple(x):i for i,x in enumerate(U)}
 for qi in qmin:
  ff=P[vals==key(qi)];assert len(ff) and modrank(ff^ff[0])==d
 dif=[[mp[tuple((a-b)%D)] for b in U] for a in U]
 exact=F((rank+1)**2,D)*F(4*D-int(beta),4*D)**d
 for b in sorted(set(map(int,w))):
  unit=F(4*D-b,4*D)**d/D
  if D*D*unit<exact:ub=D
  elif rank==4 and b==150:
   z=next(i for i,v in enumerate(w)if v>b);assert len({tuple((t*U[z])%D)for t in range(5)})==5;ub=50
  else:
   ne=[i for i in range(1,D)if w[i]<=b];g=[sum((1<<j) for j,t in enumerate(ne) if c!=t and w[dif[c][t]]<=b)for c in ne];ub=1+maxclq(g)[0]
  assert ub**2*unit<=exact
  print('carrier',d,b,ub)
