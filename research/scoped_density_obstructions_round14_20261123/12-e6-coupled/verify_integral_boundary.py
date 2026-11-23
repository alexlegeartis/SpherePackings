# Run with argv[1]=3, in the same directory as verify_cone.py.
from verify_cone import *
from sympy import Matrix, linsolve
import math
cases_checked=0
for T,(d,X) in masks.items():
  if X!=(F(3,2),F(11,7),F(10,7)) or F(64,9)*d**3>=45:continue
  S=list(T)
  if cap(S)<45:continue
  cases_checked+=1
  s=sz[S]; r=-np.unique(M[1:,S],axis=0)
  A=np.r_[r,np.diag(s)]; b=np.r_[np.full(len(r),2),s]
  Y=[(1,F(3419,4308)),(7,F(1159,4308)),(10,F(4975,2154)),(16,F(5)),(45,F(4651,4308)),(52,F(1013,1436)),(57,F(2)),(81,F(1725,1436)),(97,F(989,718)),(128,F(10431,1436))]
  deficit=[F(int(s[k]))-sum(y*int(A[i,k]) for i,y in Y) for k in range(len(s))]
  assert all(z<=0 for z in deficit)
  assert 1+sum(int(b[i])*y for i,y in Y)==45
  l=math.lcm(*map(int,s))
  B0=(A*(l//s)).astype(np.int64)
  rows=[list(map(int,B0[i])) for i,y in Y]+[[1]*len(s)]
  rhs=[F(45*l*int(b[i]),2) for i,y in Y]+[F(990)]
  for i,e in enumerate(deficit):
      if e<0:
          v=[0]*len(s);v[i]=1;rows.append(v);rhs.append(F(0))
  # Active dual eliminates all strict-deficit differences. The verified
  # affine equality forces complete translation by the class g of orbit13.
  Z=list(linsolve((Matrix(rows),Matrix(rhs))))[0]
  k=13
  assert Z[k]==45 and len(P[S[k]])==2
  g=V[P[S[k]][0]].astype(np.int64)
  ids=[i for i,e in enumerate(deficit) if e>=0]
  pool={0}
  for i in ids:pool.update(map(int,P[S[i]]))
  def quotient(x):
      x=np.asarray(x,dtype=np.int64)
      return min(int(((x+j*g)%3)@pw) for j in range(3))
  nbr={quotient(V[p]) for p in pool if all(int(((V[p]+j*g)%3)@pw) in pool for j in range(3))}
  assert 0 in nbr
  adj=[[] for _ in sorted(nbr)[1:]]
  Nl=sorted(nbr)[1:]
  lookup=set(Nl)
  for i,a in enumerate(Nl):
      for j,b in enumerate(Nl[:i]):
          if quotient(V[a]-V[b]) in lookup:
              adj[i].append(j);adj[j].append(i)
  B=[sum(1<<k for k in e) for e in adj]
  print("quotient graph",len(B),sum(map(len,adj))//2,"need clique",14,flush=True)
  # Greedy proper coloring of each induced subgraph + enumeration.
  nodes=0
  def bits(x):
      while x:
          p=x&-x;x^=p
          yield p.bit_length()-1
  def search(P,k):
      global nodes
      nodes+=1
      if k==0:return True
      if P.bit_count()<k:return False
      W=P;order=[];bounds=[];c=0
      while W:
          c+=1;A=W
          while A:
              v=A&-A; i=v.bit_length()-1
              W^=v;A^=v;A&=~B[i]
              order.append(i);bounds.append(c)
      for v,u in reversed(list(zip(order,bounds))):
          if u<k:return False
          if search(P&B[v],k-1):return True
          P^=1<<v
      return False
  assert not search((1<<len(B))-1,14)
  print("quotient excludes 45 sites, exact coloring nodes",nodes,flush=True)

assert cases_checked==1
