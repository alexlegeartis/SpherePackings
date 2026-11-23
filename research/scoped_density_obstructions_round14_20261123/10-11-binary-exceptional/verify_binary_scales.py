import numpy as np, math,itertools as I
from scipy.optimize import linprog
from fractions import Fraction as F
from collections import defaultdict
from time import time 
def diag(t):
 k,n=t[0],int(t[1:]); edges=[(i,i+1) for i in range(n-1)]
 if k=='D':edges[-1]=(n-3,n-1)
 if k=='E': edges=[(i,i+1) for i in range(n-2)]+[(2,n-1)]
 G=np.eye(n,dtype=int)*2
 for i,j in edges:G[i,j]=G[j,i]=-1
 det=n+1 if k=='A' else 4 if k=='D' else 9-n
 return n,G,det
def data(t,q):
 n,G,det=diag(t); N=q**n; powq=[q**i for i in range(n)]
 R=np.linalg.cholesky(G).T
 z=np.zeros(n,dtype=int); roots=[]
 def rec(i,r):
  if i<0:
   q0=int(z@G@z)
   if 0<q0<=8: roots.append(z.copy())
   return
  a=R[i,i+1:]@z[i+1:]; b=R[i,i]; sr=math.sqrt(max(0,r))
  for j in range(math.ceil((-sr-a)/b-1e-8),math.floor((sr-a)/b+1e-8)+1):
   z[i]=j;rec(i-1,r-(a+b*j)**2)
 rec(n-1,8)
 

 vs=np.asarray([[(i//h)%q for h in powq] for i in range(N)],dtype=int)
 def orb(dual):
  todo=set(range(N)); out=[]
  maps=[]
  for i in range(n):
   if dual:a=vs-vs[:,i,None]*G[i]
   else:
    a=vs.copy()
    a[:,i]-=vs@G[:,i]
   maps.append((a%q)@powq)
  while todo:
   g=todo.pop(); S=[g]
   for s in S:
    for a in maps:
     v=int(a[s])
     if v in todo: todo.remove(v);S.append(v)
   out.append(S)
  return out
 p=orb(False); Q=orb(True)
 V=vs[[x[0] for x in p]].copy()
 V[V>q/2]-=q
 R=np.asarray(roots); U=R@G; a=q*np.sum(U*R,axis=1)
 while 1:
  gain=2*V@U.T-a
  im=np.argmax(gain,axis=1); b=gain[np.arange(len(V)),im]>0
  if not np.any(b):break
  V[b]-=q*R[im[b]]
 assert np.all((V%q)@powq==[x[0] for x in p])
 depths=np.sum((V@G)*V,axis=1)

 lut={2:[2,-2],3:[2,-1,-1],4:[2,0,-2,0],6:[2,1,-1,-2,-1,1]}[q]; lut=np.asarray(lut,dtype=int)
 M=np.array([[lut[(vs[y]@vs[S].T)%q].sum() for S in p] for y in [T[0] for T in Q]])
 assert np.max(np.abs(M-np.rint(M)))<1e-6
 return n,det,depths,np.rint(M).astype(int),np.array([len(S) for S in p],int)

T=['A'+str(i) for i in range(1,10)]+['D'+str(i) for i in range(4,10)]+['E6','E7','E8']
data0={}
def kron(a):
 z=a[0]
 for x in a[1:]:z=np.kron(z,x)
 return z
def det(m):
 if len(m)==1:return int(m[0][0])
 return sum((-1)**j*int(m[0][j])*det([r[:j]+r[j+1:] for r in m[1:]]) for j in range(len(m)))
for n in [10,11]:
 b2={10:F(1,768),11:F(1,1024)}[n];bs2={10:F(5,128)**2,11:F(9,256)**2}[n]
 for r in [2]:
  for ss in I.combinations_with_replacement(T,r):
   if sum(int(s[1:]) for s in ss)!=n or 'E8' in ss:continue
   if ss not in [('D4','E6'),('D4','E7')]:continue
   for qs in [(2,2)]:
    for s,q in zip(ss,qs):
     if (s,q) not in data0:data0[(s,q)]=data(s,q)
    ds=[data0[(s,q)] for s,q in zip(ss,qs)]
    D=np.array(list(I.product(*[l[2] for l in ds])),int)
    M=kron([l[3] for l in ds]); S=kron([l[4] for l in ds])
    den0=math.lcm(*(q*q for q in qs)); D*=np.array([den0//(q*q) for q in qs]); lim=2*den0
    B=np.r_[D[1:],np.eye(r,dtype=int)*lim].tolist()
    _,idx=np.unique(B,axis=0,return_index=True)
    B=[B[i] for i in sorted(idx)]
    B=[b for b in B if sum(b)<=lim or max(b)==lim and b.count(0)==r-1]
    V=set()
    for b in I.combinations(B,r):
     h=det(b)
     if h==0:continue
     v=tuple(F(det([row[:i]+[lim]+row[i+1:] for row in b]),h) for i in range(r))
     if min(v)>=1:V.add(v)
    polys={}
    for v in V:
     vol2=math.prod(l[1]*x**l[0] for l,x in zip(ds,v)); target=math.sqrt(float(b2*vol2)*2**n)
     den=math.lcm(*(x.denominator for x in v))
     mask=tuple(np.flatnonzero(D@np.array([int(x*den) for x in v])>=lim*den))
     if F(sum(map(int,S[list(mask)]))+1)**2 < b2 * vol2 * 2**n:continue
     if mask not in polys or polys[mask][0]>vol2:polys[mask]=vol2,v,target
    res=[]
    for mask,(vol2,v,target) in polys.items():
     R=-M[1:,mask]
     R=np.unique(R,axis=0); A=np.r_[R,np.diag(S[list(mask)])]; YB=np.r_[np.full(len(R),2**r),S[list(mask)]]
     z=linprog(-S[list(mask)],A_ub=A,b_ub=YB,bounds=(0,None), method='highs')
     assert z.success
     Y=[F(max(0,-x)).limit_denominator(10**6) for x in z.ineqlin.marginals]
     sparse=[(j,y) for j,y in enumerate(Y) if y]
     corr=sum(max(F(0),F(int(s))-sum(y*int(A[j,i]) for j,y in sparse)) for i,s in enumerate(S[list(mask)]))
     cap=1+sum(int(b)*y for b,y in zip(YB,Y))+corr
     def floor(F0):return F0.numerator//F0.denominator
     ratio=F(floor(cap)**2)/(b2*vol2*2**n)
     assert ratio<=1, (ss,qs,v,cap,ratio)
     if ratio>=1: print('AT-LEAST-LATTICE',ss,qs,v,cap,ratio,flush=True)
     res.append((ratio,v,cap,target,len(mask)))
    print(ss,qs,len(V),len(polys),[(float(z[0])**.5,tuple(map(str,z[1])),str(z[2]),z[3],z[4]) for z in sorted(res,reverse=True)[:3]],flush=True)
