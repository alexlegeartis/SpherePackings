from even import *
from anchored_replay import N
import json
from math import prod
from anchored_replay import COEF,q, SZ as Z10
S12={2:1331791267165632000,6:618651644974272000,10:123081209203680000,14:8636968424448000,18:143169114816000,22:231390432000}
S12|={-k:v for k,v in S12.items()}
N12=sum(S12.values())
def valid(p,c,s):
 m=(14+s)//2; n=(42+s)//2
 G=((10,0,p,c),(0,10,c,p+c),(p,c,14,m),(c,p+c,m,n))
 a=140-p*p-c*c;b=10*n-c*c-(p+c)**2;d=10*m-p*c-c*(p+c)
 return a>=0 and b>=0 and a*b>=d*d and not any(0<q(G,t)<10 for t in COEF)
def images(p,c,s):
 out=[(-p,-c,s)]
 if s==-14:out.append((c,p+c,14))
 if s==14:out.append((c-p,p,-14))
 for u,v in product(range(-1,2),repeat=2):
  if not(u or v) or 5*(u*u+v*v)+p*u+c*v: continue
  h=-10*u*u+40*u*v+10*v*v+2*((2*c-p)*u+(2*p+c)*v)
  out.append((p+10*u,c+10*v,s+h))
 return out
def orbits():
 A={(p,c,s)for p,c,s in product(range(-7,8),range(-7,8),S12) if valid(p,c,s)}
 B=[];forbid=[]
 while A:
  x=min(A);todo=[x];Z={x}
  while todo:
   for y in images(*todo.pop()):
    if y not in Z:todo.append(y);Z.add(y)
  (B if Z<=A else forbid).append(tuple(sorted(Z&A)))
  A-=Z
 return B,forbid

def shellrows(B):
 R=[];V=[];M=[]
 def add(f,t,l):
  R.append([sum(f(*a) for a in bb) for bb in B]);V.append(F(t));M.append(l)
 for s,A in S12.items():
  add(lambda p,c,h:int(h==s), A,('size',s))
  add(lambda p,c,h:(5*p*p+(2*c-p)**2)*int(h==s),F(700-10*s,48)*A,('cov1',s))
  add(lambda p,c,h:p*(2*c-p)*int(h==s),F(5*s-70,48)*A,('cov2',s))
 for i in range(11):
  for j in range(11-i):
   if (i+j)%2 or not (i+j):continue
   t=0
   if not i%2 and not j%2:
    m=(i+j)//2
    t=F(N12*140**m*prod(range(1,i,2))*prod(range(1,j,2)),prod(96+2*h for h in range(m)))
   add(lambda p,c,s:p**i*c**j,t,('mom',i,j))
 return R,V,M
