"""Standalone replay; Python 3, NumPy, and SciPy. See README.md."""

import os
for s in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS']: os.environ[s]='1'
import sys,itertools as it,math,time,json
from collections import Counter,defaultdict
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
def pr(*a): print(*a,flush=True)
def gdet(a,b,c):return a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0])
def cross(a,b):return tuple(a[(i+1)%3]*b[(i+2)%3]-a[(i+2)%3]*b[(i+1)%3] for i in range(3))
def vclasses(k):
 pts=[x+((-sum(x))%3,) for x in it.product(range(3),repeat=k)]
 def key(x):return min(tuple(sorted(x)),tuple(sorted(-u%3 for u in x)))
 groups=defaultdict(list)
 for x in pts:groups[key(x)].append(x)
 groups=sorted(groups.values(),key=lambda g:(any(g[0]),key(g[0])))
 labels=[]; rr=[]; rows=[]
 for g in groups:
  x=g[0]; vals=[sum((3*t+a)**2 for a,t in zip(x,b)) for b in it.product([-2,-1,0,1],repeat=k+1) if sum(x)+3*sum(b)==0]; labels.append(2*min(vals))
 for u in it.product(range(3),repeat=k):
  v=tuple(sum(2 if sum(p*q for p,q in zip(u,x))%3==0 else -1 for x in g)//2 for g in groups);rr.append(v)
 typ=sorted(Counter(rr).items(),key=lambda t:sum(t[0]),reverse=True)
 for j,g in enumerate(groups):
  row=[m*e[j]//len(g) for e,m in typ];assert all(row[i]*len(g)==m*e[j] for i,(e,m) in enumerate(typ));rows.append(row)
 return labels,np.array(rows,dtype=np.int64)
KC=[0,12,18,24];K=np.array([[1,252,224,252],[1,9,8,-18],[1,9,-19,9],[1,-18,8,9]],dtype=np.int64)
def scan(k):
 vc,V=vclasses(k);C=list(it.product(KC,KC,vc));O=np.kron(np.kron(K,K),V); mul=O[0,1:]; CM=O[:,1:]/mul; group=3**(12+k)
 assert tuple(O[0])==tuple(np.concatenate(([1],mul))) and sum(O[0])==group
 assert all(sum(r)==0 for r in O[1:]); D=10**14; memo={}
 def lp(inds):
  key=tuple(inds)
  if key not in memo:
   A=CM[list(key),:]; res=linprog(np.ones(len(mul)),A_ub=A,b_ub=-np.ones(len(key)),method='highs')
   if not res.success:raise ValueError(res.message)
   s0=[max(0,round(D*x/int(m))) for x,m in zip(res.x,mul)]; rc=max([0]+[D+sum(int(b)*c for b,c in zip(O[j,1:],s0)) for j in key]); val=D+sum(int(b)*c for b,c in zip(mul,s0))+rc*(group-1)
   assert min(s0)>=0; memo[key]=(val//D,val/D)
  return memo[key][0]
 rows=sorted(set(C[1:])|{(36,0,0),(0,36,0),(0,0,36)})
 scenarios={}; nv=0
 for a,b,c in it.combinations(rows,3):
  d=gdet(a,b,c)
  if not d:continue
  u,v,w=cross(a,b),cross(b,c),cross(c,a); num=[36*(u[j]+v[j]+w[j]) for j in range(3)]
  if d<0:d=-d;num=[-n for n in num]
  if min(num)<d or num[0]>num[1]:continue
  common=math.gcd(d,math.gcd(*num)); ts=tuple(x//common for x in num)+(d//common,);scenarios[ts]=None
 pr('start',k,'metrics',len(scenarios),'vc',vc)
 passing=[]; lpassing=[]; top=[]; counts=Counter()
 benchmark2=[0,F(1,2),F(1,3),F(1,2)][k]; lat2=[0,F(1,2),F(1,3),F(1,3)][k]; 
 for (a,b,c,d) in scenarios:
  ids=tuple(j for j,(i1,i2,iv) in enumerate(C[1:],1) if a*i1+b*i2+c*iv>=36*d)
  M=lp(ids)
  den2=729**2*(F(a*b,d*d))**12*(2**k)*(k+1)*F(c,d)**k
  opt2=F(M*M,1)/den2
  p3=3**max(j for j in range(13+k) if 3**j<=M)
  q2=F(p3*p3,1)/den2; tr=tuple(str(F(x,d)) for x in (a,b,c)); entry=(float(opt2)**.5,tr,M,p3,float(q2)**.5)
  top.append(entry)
  if opt2>benchmark2:passing.append(entry)
  if q2>lat2:lpassing.append(entry)
 pr('result',k,'LPs',len(memo),'nonlinear left',len(passing),'linear left',len(lpassing),'max',sorted(top,reverse=True)[:7]);pr('lin',k,sorted(lpassing,reverse=True)[:40]);pr('nonlin',k,sorted(passing,reverse=True)[:24]);return passing,lpassing

results={k:scan(k) for k in (1,2,3)}
v=np.array(list(it.product(range(3),repeat=6)),dtype=np.int16); key=np.array([3**i for i in range(5,-1,-1)],dtype=np.int16)
q=(v*v).sum(1)%3;wws=[np.where(q==1,4,np.where(q==2,8,6)),np.where(q==1,8,np.where(q==2,4,6))];wws[0][0]=wws[1][0]=0
v3=np.array(list(it.product(range(3),repeat=3)),dtype=np.int16); powers3=np.array([9,3,1]); GL=[]
def rankF(M):
 a=np.array(M,dtype=int)%3;h=0
 for j in range(a.shape[1]):
  nz=np.flatnonzero(a[h:,j])
  if not len(nz):continue
  i=h+nz[0];a[[i,h]]=a[[h,i]];a[h]=a[h]*a[h,j]%3;ix=np.flatnonzero(a[:,j]);ix=ix[ix!=h];a[ix]=(a[ix]-a[ix,j,None]*a[h])%3;h+=1
  if h==len(a):break
 return h
GL=np.array([np.array([a,b,c]) for a in v3[1:] for b in v3[1:] for c in v3[1:] if gdet(a,b,c)%3],dtype=np.int16)
perms3=((v3@GL)%3)@powers3
def formVals(A):return ((v3@A%3)**2).sum(axis=1)%3
def embed(Q):
 G=np.array([[(Q[((np.eye(3,dtype=int)[i]+np.eye(3,dtype=int)[j])%3)@powers3]-Q[3**(2-i)]-Q[3**(2-j)])*2%3 for j in range(3)] for i in range(3)])
 def rec(rows):
  i=len(rows)
  if i==3:return np.array(rows)
  S=(q==G[i,i]);nz=1
  for j,r in enumerate(rows):S&= ((v@r)%3==G[i,j])
  for r in v[S]:
   if rankF(rows+[r])!=i+1:continue
   w=rec(rows+[r])
   if w is not None:return w
  return None
 return rec([])
def formpairs():
 E=np.stack([v3[:,i]*v3[:,j] for i,j in it.combinations_with_replacement(range(3),2)],axis=1)
 forms=np.array([(E@np.array(a))%3 for a in it.product(range(3),repeat=6)],dtype=np.int16)
 done=set(); pairs=[]
 for qa in forms:
  t=qa.tobytes()
  if t in done:continue
  shifts=qa[perms3]; done.update(x.tobytes() for x in shifts);A=embed(qa)
  if A is None:continue
  auto=perms3[np.all(shifts==qa,axis=1)];doneb=set()
  for qb in forms[np.all(((forms!=1)|(qa[None,:]==2))&((qa[None,:]!=1)|(forms==2)),axis=1)]:
   tb=qb.tobytes()
   if tb in doneb:continue
   doneb.update(x.tobytes() for x in qb[auto]);B=embed(qb)
   if B is None:continue
   assert all(formVals(A)==qa) and all(formVals(B)==qb); pairs.append((A,B,qa,qb))
 pr('rank3 forms',len(pairs),sum(np.all((a+b)%3==0) for A,B,a,b in pairs));return pairs
pairs3=formpairs()
def ds3(A,B):
 R=A.tolist();extra=[]
 for e in np.eye(6,dtype=int):
  if rankF(R+[e])>len(R):R.append(e);extra.append(e)
 X=v3@A%3; Y=v3@B%3; Z=v3@np.array(extra)%3
 ds=np.full((27,729),100,dtype=int);bad=np.zeros_like(ds,dtype=bool)
 for a,b in zip(X,Y):
  ai=((a+Z)%3)@key;bi=((b+v)%3)@key; ds=np.minimum(ds,wws[0][ai,None]+wws[0][None,bi]);bad|=(q[ai,None]+q[None,bi])%3!=0
 return ds.ravel(),bad.ravel()
data3=[ds3(A,B) for A,B,qa,qb in pairs3];pr('neighbors',[(int(sum(d>=12)),int(sum((d>=12)&bad))) for d,bad in data3])
v2=np.array(list(it.product(range(3),repeat=2)),dtype=np.int16); powers2=np.array([3,1]); GL=[]
def rankF(M):
 a=np.array(M,dtype=int)%3;h=0
 for j in range(a.shape[1]):
  nz=np.flatnonzero(a[h:,j])
  if not len(nz):continue
  i=h+nz[0];a[[i,h]]=a[[h,i]];a[h]=a[h]*a[h,j]%3;ix=np.flatnonzero(a[:,j]);ix=ix[ix!=h];a[ix]=(a[ix]-a[ix,j,None]*a[h])%3;h+=1
  if h==len(a):break
 return h
GL=np.array([np.array([a,b]) for a in v2[1:] for b in v2[1:] if (a[0]*b[1]-a[1]*b[0])%3],dtype=np.int16)
perms2=((v2@GL)%3)@powers2
def formVals2(A):return ((v2@A%3)**2).sum(axis=1)%3
def embed2(Q):
 G=np.array([[(Q[((np.eye(2,dtype=int)[i]+np.eye(2,dtype=int)[j])%3)@powers2]-Q[3**(1-i)]-Q[3**(1-j)])*2%3 for j in range(2)] for i in range(2)])
 def rec(rows):
  i=len(rows)
  if i==2:return np.array(rows)
  S=(q==G[i,i]);nz=1
  for j,r in enumerate(rows):S&= ((v@r)%3==G[i,j])
  for r in v[S]:
   if rankF(rows+[r])!=i+1:continue
   w=rec(rows+[r])
   if w is not None:return w
  return None
 return rec([])
def formpairs2():
 E=np.stack([v2[:,i]*v2[:,j] for i,j in it.combinations_with_replacement(range(2),2)],axis=1)
 forms=np.array([(E@np.array(a))%3 for a in it.product(range(3),repeat=3)],dtype=np.int16)
 done=set(); pairs=[]
 for qa in forms:
  t=qa.tobytes()
  if t in done:continue
  shifts=qa[perms2]; done.update(x.tobytes() for x in shifts);A=embed2(qa)
  if A is None:continue
  auto=perms2[np.all(shifts==qa,axis=1)];doneb=set()
  for qb in forms[np.all(((forms!=1)|(qa[None,:]==2))&((qa[None,:]!=1)|(forms==2)),axis=1)]:
   tb=qb.tobytes()
   if tb in doneb:continue
   doneb.update(x.tobytes() for x in qb[auto]);B=embed2(qb)
   if B is None:continue
   assert all(formVals2(A)==qa) and all(formVals2(B)==qb); pairs.append((A,B,qa,qb))
 pr('rank2 forms',len(pairs),sum(np.all((a+b)%3==0) for A,B,a,b in pairs));return pairs
pairs2=formpairs2()
def ds2(A,B):
 R=A.tolist();extra=[]
 for e in np.eye(6,dtype=int):
  if rankF(R+[e])>len(R):R.append(e);extra.append(e)
 X=v2@A%3; Y=v2@B%3; Z=np.array(list(it.product(range(3),repeat=4)),dtype=np.int16)@np.array(extra)%3
 ds=np.full((81,729),100,dtype=int);bad=np.zeros_like(ds,dtype=bool)
 for a,b in zip(X,Y):
  ai=((a+Z)%3)@key;bi=((b+v)%3)@key; ds=np.minimum(ds,wws[0][ai,None]+wws[0][None,bi]);bad|=(q[ai,None]+q[None,bi])%3!=0
 return ds.ravel(),bad.ravel()
data2=[ds2(A,B) for A,B,qa,qb in pairs2];pr('neighbors',[(int(sum(d>=12)),int(sum((d>=12)&bad))) for d,bad in data2])
def charts(A):
 R=A.tolist();extra=[]
 for e in np.eye(6,dtype=int):
  if rankF(R+[e])>len(R):R.append(e);extra.append(e)
 Z=np.array(list(it.product(range(3),repeat=4)),dtype=np.int16)@np.array(extra)%3
 X=v2@A%3;ls=((Z[:,None,:]+X[None,:,:])%3)@key; addr=np.empty(729,int);addr[ls.ravel()]=np.arange(729);return Z,X,addr
def revmap(A,B):
 Z,X,ma=charts(A);T,Y,mb=charts(B);return ((mb//9)[None,:]*729+((Z[:,None,:]-X[mb%9][None,:,:])%3)@key).ravel()
def nfft(x,inverse=False):
 p=1009; w=next(i for i in range(2,p) if i**3%p==1)
 if inverse:w=(w*w)%p
 b=x.reshape((len(x),)+(3,)*6).astype(np.int64)
 for j in range(1,7):
  h=b.swapaxes(1,j);z=np.stack([(h[:,0]+h[:,1]+h[:,2])%p,(h[:,0]+w*h[:,1]+(w*w%p)*h[:,2])%p,(h[:,0]+(w*w%p)*h[:,1]+w*h[:,2])%p],axis=1);b=z.swapaxes(1,j)
 if inverse:b=b*pow(729,-1,p)%p
 return b.reshape(len(x),729)
def exactrows(ds):
 x=(ds.reshape(81,729)>=12).astype(int);f=nfft(x);g=nfft(x,True)*729%1009;rs=np.array(list(it.product(range(3),repeat=4)));kk=np.array([27,9,3,1]);out=np.zeros_like(x)
 for j,r in enumerate(rs):
  J=((rs-r)%3)@kk;co=nfft(f*g[J]%1009,True);assert int(np.max(co))<=729;out[j]=(co>0).sum(0)
 return out.ravel()

rem2={}
for i,((A,B,qa,qb),(ds,bad)) in enumerate(zip(pairs2,data2)):
 rr=revmap(A,B);dst,_=ds2(B,A);assert np.all(ds==dst[rr])
 cap=np.minimum(exactrows(ds),exactrows(dst)[rr])
 allow=(ds>=12)&((~bad)|(cap>=53))
 if np.any(allow&bad):rem2[i]=allow
 pr("hard defect count",i,sum(allow&bad))
v5=np.array(list(it.product(range(3),repeat=5)),np.int16);sub5=((v5[:,None,:]-v5[None,:,:])%3)@(3**np.arange(4,-1,-1));addr=np.arange(59049)
def diffone(s):return sub5[addr//243,s//243]*243+sub5[addr%243,s%243]
def clq(L,s,goal=53):
 ix=np.flatnonzero(L&L[diffone(s)]);h=ix[:,None]//243;l=ix[:,None]%243
 D=L[243*sub5[h,h.T]+sub5[l,l.T]];A=[int.from_bytes(row.tobytes(),'little') for row in np.packbits(D,axis=1,bitorder='little')]
 calls=0
 def visit(P,n):
  nonlocal calls
  calls+=1
  if n==0:return True
  # first-fit disjoint stable greedy layers, retain order for monotone bounds
  order=[];left=P;c=0
  while left:
   c+=1;rem=left
   while rem:
    z=rem&-rem;j=z.bit_length()-1;order.append((z,j,c));left^=z;rem=(rem^z)&~A[j]
  for z,j,c in reversed(order):
   if c<n:return False
   if visit(P&A[j],n-1):return True
   P^=z
  return False
 found=visit((1<<len(ix))-1,goal);return found,calls
def allbb():
 for i,L in rem2.items():
  bad=data2[i][1];cc=0;mx=0;k=0
  for ss in np.flatnonzero(L&bad):
   b,n=clq(L,int(ss));assert not b;cc+=n;mx=max(mx,n);k+=1
  pr('all defect cliques',i,k,cc,mx)
allbb()

marks1={0:[1,1,1,0,0,0],1:[1,0,0,0,0,0],2:[1,1,0,0,0,0]}
T=np.zeros((6,6),np.int16)
for j in range(0,6,2):T[j:j+2,j:j+2]=[[1,1],[1,-1]]
full=(v@T%3)@key
for m in [marks1[0],marks1[1],marks1[2]]:
 m=np.array(m); a=q[m@key];w=wws[0]+wws[0][((v@T+m)%3)@key]
 pr('one point trades',a,dict(Counter(map(int,w))),sum(w>=12)+1)
def localbad(ma):
 m=np.array(marks1[ma]);X=(v@T)%3;h=(X@m)%3;J=(wws[0]+wws[0][((X+m)%3)@key])<12;out={}
 for k in range(3):
  mask=J&(h==k);dd=sum(mask);pd=v@v[mask].T%3; ar=np.where(pd==0,2,-1).sum(1);ai=np.where(pd==1,1,np.where(pd==2,-1,0)).sum(1);sp=(ar**2+3*ai**2)//4
  # remove all 3 dual characters that are constant on the hyperplane
  local=(v@v[h==0].T)%3; ss=int(np.max(sp[np.any(local!=0,axis=1)]));out[k]=(int(dd),ss)
 pr('slices',ma,out);return out
slicevals=[localbad(i) for i in range(3)]
def collisions():
 # use labels in first coordinates: an incompatible point is (a,T(a+s)): forbidden first phases a+E_s, with E_s from w(x)+w(T(x-s))<4.
 W=wws[0];Wt=wws[0][(v@T%3)@key];E=np.array([(W+Wt[((v-z)%3)@key]<12) for z in v]); f=nfft(E.astype(int));gf=nfft(E.astype(int),True)*729%1009
 for mark in marks1.values():
  h=int(np.array(mark)@key);corr=nfft(f[h,None,:]*gf%1009,True); ok=np.ones_like(corr,dtype=bool);ok[0]=False
  # at offset a, intersection E_h against a+E_s; permitted if (a,T(a+s-h)) satisfies intra-layer minimum.
  Aplus=((v[None,:,:]+v[:,None,:]-v[h])%3)@key;ok&=(W[None,:]+Wt[Aplus])>=12
  axis=np.array([j not in (h,(2*v[h]%3)@key) for j in range(729)]);pr('by axis',q[h],int(np.max(corr[ok&axis[:,None]])),Counter(map(int,np.max(np.where(ok,corr,-1),axis=1))));pr('all two outsiders',int(q[h]),int(np.max(corr[ok])),int(np.min((E.sum(1)[:,None]+E[h].sum()-corr)[ok])))
collisions()
for ma in range(3):
 m=np.array(marks1[ma]);X=(v@T)%3;J=(wws[0]+wws[0][((X+m)%3)@key])<12;F=(wws[0]+wws[0][((X-m)%3)@key])<12
 cc=[int(sum(J&F[((v+b)%3)@key])) for b in v];pr('opposing overlap',ma,max(cc))
 for mm in (1,2,3,4,8,12,20,24):
  vall=[(d,s) for d,s in slicevals[ma].values() if d>2]
  B=243+sum(min(243-d,(243*s*(243-mm))//(d*d*mm+s*(243-mm))) for d,s in vall)+6*mm;pr(ma,mm,B)

assert max(243+sum(min(243-dd,(243*ss*(243-x))//(dd*dd*x+ss*(243-x))) for dd,ss in slicevals[ma].values() if dd>2)+6*x for x in range(2,25) for ma in range(3))==470
