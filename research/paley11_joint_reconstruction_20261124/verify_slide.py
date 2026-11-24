from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
def read(name):
 return [tuple(map(int,s.split())) for s in Path(name).read_text().splitlines()]
def verify(name,K=1,B=2):
 X=read(name);P=4*K
 def canon(x):
  a=[u%P for u in x]
  for b in range(B):a[4*b]+=P*(sum(u//P for u in x[4*b:4*b+4])%2)
  return tuple(a)
 R=set(map(canon,X));assert len(R)==len(X)
 def metric(x,y):
  d=[a-b for a,b in zip(x,y)];v=[(u+P//2)%P-P//2 for u in d]
  D=sum(t*t for t in v)
  for b in range(B):
   sl=slice(4*b,4*b+4)
   d0=sum((x-y)//P for x,y in zip(d[sl],v[sl]))%2
   D+=d0*(P*P-2*P*max(map(abs,v[sl])))
  return D
 A=[metric(x,y)for x,y in combinations(X,2)]
 if name.startswith("slide20"):
  M={r:a for r,a in read("slide20-axes.txt")}
  assert len(M)==4 and all(8<=a<11 for a in M.values())
  for i,j in combinations(range(len(X)),2):
   ax={M[r] for r in (i,j) if r in M}
   if ax:
    u,v=X[i],X[j]
    dropped=sum(((u[a]-v[a]+P//2)%P-P//2)**2 for a in ax)
    assert metric(u,v)-dropped>=16*K*K
 assert min(A)==16*K*K
 anchor=X[0];stab=0
 for x in X:
  t=tuple(a-b for a,b in zip(x,anchor))
  if {canon(tuple(a+b for a,b in zip(u,t)))for u in X}==R:stab+=1
 frac=sum(d%(4*K*K)!=0 for d in A)
 vol=2**(11+B)
 assert F(len(X),vol)==F(9,256)
 # Pure axis pair contacts (or same-site period) imply a first-order trace obstruction even with joint site motion.
 for j in range(11):
  assert any(canon(u[:j]+(u[j]+P,)+u[j+1:]) in R for u in X)
 return stab,frac
if __name__=='__main__':
 S=read('trial-1-288.txt');assert all(all(a%2==0 for a in u[8:])for u in S)
 assert all(all(len({a%2 for a in u[b:b+4]})==1 for b in (0,4))for u in S)
 assert len(set(read("atlas-tie1.txt"))-set(S))==24
 T=[list(x) for x in read("slide20-288-k257.txt")]
 for j,(r,a) in enumerate(read("slide20-axes.txt"),1): T[r][a]-=2*j
 assert all(x%257==0 for u in T for x in u)
 assert len({tuple(x//257 for x in u) for u in T}-set(S))==20
 for name,K,B,target in [('atlas-tie1.txt',1,2,1),('slide20-288-k257.txt',257,2,1),('atlas-tie0.txt',1,1,1)]:
  s,f=verify(name,K,B);assert s==target
  if B==2:assert f>0
  print(name,'primitive',s,'fractional distance pair count',f)
