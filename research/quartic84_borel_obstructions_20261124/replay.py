from fractions import Fraction
import json,itertools
from verify import row,mul,transpose,q,det,pow,identity,A,G,R,J
from receipt import read
p=43;n=84
def rref(M,p):
 c=[list(w)for w in M];i=0
 for j in range(len(c[0])):
  a=next((a for a in range(i,len(c))if c[a][j]%p),None)
  if a is None:continue
  c[i],c[a]=c[a],c[i]
  r=pow0(c[i][j]%p,-1,p);c[i]=[x*r%p for x in c[i]]
  for k,w in enumerate(c):
   if k==i:continue
   x=w[j]%p
   if x:c[k]=[(t-x*v)%p for t,v in zip(w,c[i])]
  i+=1
  if i==len(c):break
 return tuple(map(tuple,c[:i]))
pow0=__builtins__['pow'] if isinstance(__builtins__,dict) else __builtins__.pow
D=read['borel'];Rt,Jt=D['R43'],D['J43'];Id=identity(n)
assert len(rref(A,p))==64
assert mul(mul(Rt,A),transpose(Rt))==A and pow(Rt,21)==Id
assert mul(Rt,Jt)==mul(Jt,Rt) and mul(Jt,Jt)==[[-x for x in r]for r in Id]
assert mul(mul(Jt,A),transpose(Jt))==A
S={}
for (l,m),(X,Y) in zip(D['eigenpairs'],D['eigenbases']):
 assert l*m%p==1 and l!=m
 for z,C in [(l,X),(m,Y)]:
  assert z not in S and pow0(z,21,p)==1
  assert all(all(x%p==0 for x in r)for r in mul(C,A))
  assert [[x%p for x in r]for r in mul(C,Rt)]==[[z*x%p for x in r]for r in C]
  assert len(rref(C,p))==2 and rref(mul(C,Jt)+C,p)==rref(C,p)
  S[z]=C
 Z0=mul(mul(X,A),transpose(Y))
 assert all(x%p==0 for r in Z0 for x in r)
 assert det([[x//p for x in r]for r in Z0])%p
assert len(S)==10 and len(rref(sum(S.values(),[]),p))==20 and pow0(p-1,(p-1)//2,p)==p-1
masks=set()
for i,x in enumerate(D['variants']):
 h=x['H'];b=x['mask'];v=x['witness']['v43'];N=x['witness']['q']
 assert tuple(b) not in masks;masks.add(tuple(b))
 c=sum((S[pa[j]]for pa,j in zip(D['eigenpairs'],b)),[])
 assert rref(h,p)==rref(c,p) and len(h)==10
 Z0=mul(mul(h,A),transpose(h))
 assert not any(x%(p*p)for r in Z0 for x in r) and not any(Z0[j][j]%(2*p*p)for j in range(10))
 cols=[next(i for i,u in enumerate(r)if u)for r in h]
 assert all((a-b)%p==0 for a,b in zip(v,row([v[i]%p for i in cols],h)))
 assert N<=10 and q(v,A)==N*p*p
B=D['canonical_numerator_basis']
assert mul(mul(B,A),transpose(B))==[[p*p*x for x in r]for r in G]
assert mul(B,Rt)==mul(R,B) and mul(B,Jt)==mul(J,B)
h=D['variants'][D['canonical_index']]['H']
assert rref(B,p)==rref(h,p) and abs(det(B))==p**74
print('PASS all',len(masks),'R21-stable Gaussian maximal isotropic Borel glues, each with a literal q<=10')
K=read['binary'];assert K['g']==7
T=pow(R,3);xs=K['witnesses'];vals={}
for label,x in xs.items():
 assert q(x,G)==10
 y=[a%2 for a in x];Z=y.copy()
 for _ in range(1,7):
  y=[a%2 for a in row(y,T)]; Z=[a^b for a,b in zip(Z,y)]
 label0=sum(Z[k]<<i for i,k in enumerate(K['cols']))
 assert label0==int(label)
 vals[int(label)]=Z
U={0:[0]*84}
for i in K['basis']:U|={k^i:[a^b for a,b in zip(z,vals[i])]for k,z in U.items()}
assert len(U)==32 and all(z==vals[k]for k,z in U.items()if k)
print('PASS order7 invariant binary quotients of rank<=4 blocked by 31 literal q10 projected residues')
