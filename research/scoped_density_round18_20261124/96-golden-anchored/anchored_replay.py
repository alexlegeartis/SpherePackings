"""Two golden shortest anchors in an extremal rank-96 parent."""
from fractions import Fraction as F
from itertools import product
from math import prod
from collections import defaultdict
N=565866362880
SZ={2:231390432000,6:49886928000,10:1655821440}
S=tuple(sorted(SZ.keys()|set(-x for x in SZ)))

def q(G,t):return sum(x*y*a for x,row in zip(t,G) for y,a in zip(t,row))
COEF=[t for t in product(range(-2,3),repeat=4) if any(t) and next(a for a in t if a)>0]
def G(p,c,s):
 m=(10+s)//2; n=(30+s)//2
 return ((10,0,p,c),(0,10,c,p+c),(p,c,10,m),(c,p+c,m,n))
def valid(p,c,s):
 g=G(p,c,s)
 # Schur complement on orthogonal equal-length anchors
 a=100-p*p-c*c;b=10*g[3][3]-c*c-(p+c)**2;d=10*g[2][3]-p*c-c*(p+c)
 if a<0 or b<0 or a*b<d*d:return False
 return not any(0<q(g,t)<10 for t in COEF)
def images(p,c,s):
 out=[(-p,-c,s)]
 if s==-10:out.append((c,p+c,10))
 if s==10:out.append((c-p,p,-10))
 for u in (-1,0,1):
  for v in (-1,0,1):
   if not (u or v) or 5*(u*u+v*v)+p*u+c*v:continue
   h=-10*u*u+40*u*v+10*v*v+2*((2*c-p)*u+(2*p+c)*v)
   out.append((p+10*u,c+10*v,s+h))
 return out

def orbits():
 A={(p,c,s) for p,c,s in product(range(-5,6),range(-5,6),S) if valid(p,c,s)}
 special={(-10,0,-10),(10,0,-10),(0,-10,10),(0,10,10)}
 A|=special; Bs=[]; forbidden=[]
 while A:
  x=min(A); todo=[x]; seen={x}
  while todo:
   for y in images(*todo.pop()):
    if y not in seen:seen.add(y);todo.append(y)
  Bs.append(tuple(sorted(seen))) if seen<=A else forbidden.append(tuple(sorted(seen&A)))
  A-=seen
 return sorted(Bs),forbidden

def matrix():
 B,forbid=orbits()
 rows=[]; vals=[]; names=[]
 def add(f,t,label):
  rows.append([sum(f(*a) for a in bb) for bb in B]);vals.append(F(t));names.append(label)
 for s in S:
  add(lambda p,c,h:int(h==s),SZ[abs(s)],('size',s))
 # Individual signed shell projected covariance and its mixed difference.
 # For x=(norm10, s=-10), <x_+,x_+>=5-sqrt(5), <x_-,x_->=5+sqrt(5).
 # z in shell s gives f+/- = (p +/-(2c-p)/sqrt(5))/2 with mixed covariance not needed.
 for s in S:
  A=SZ[abs(s)]
  add(lambda p,c,h:(5*p*p+(2*c-p)**2)*int(h==s),F(500-10*s,48)*A,('shell-Q',s))
  add(lambda p,c,h:p*(2*c-p)*int(h==s),F(5*s-50,48)*A,('shell-B',s))
 # Ordinary shortest shell, an 11-design in R^96 at norm 10.
 for i in range(11):
  for j in range(11-i):
   if (i+j)%2 or not (i+j):continue
   t=0
   if not i%2 and not j%2:
    m=(i+j)//2
    t=F(N*100**m*prod(range(1,i,2))*prod(range(1,j,2)),prod(96+2*h for h in range(m)))
   add(lambda p,c,s:p**i*c**j,t,('mom',i,j))
 # For odd/odd part of weight48,index(5,5), only two weak forms can
 # have order>=5. At q^5 their leading Laurent coefficient vectors
 # factor as phi_{-1,2,0}*(u-2)^a(u+10)^(3-a), a=2 or 3.
 from collections import defaultdict
 def mul(x,y):
  q=defaultdict(int)
  for a,u in x.items():
   for b,v in y.items():q[a+b]+=u*v
  return dict(q)
 def e(a):
  Q={1:1,-1:-1}
  for t in [-2]*a+[10]*(3-a):Q=mul(Q,{-1:1,0:t,1:1})
  return Q
 E={k:e(k) for k in (2,3)}
 for p,c in product(range(1,6),repeat=2):
  B1=E[2].get(p,0)*E[3].get(c,0); B2=E[3].get(p,0)*E[2].get(c,0)
  # (4,3),(3,4) give anchor basis, see computed determinant.
  a,b=E[2][4]*E[3][3],E[3][4]*E[2][3]
  cc,dd=E[2][4]*E[3][4],E[3][4]*E[2][4]
  d=a*dd-b*cc
  assert d
  A1=B1*dd-B2*cc; A2=B2*a-B1*b
  def oo(u,v,h):
   return (int((u,v)==(p,c))-int((u,v)==(p,-c)))*d-(int((u,v)==(4,3))-int((u,v)==(4,-3)))*A1-(int((u,v)==(4,4))-int((u,v)==(4,-4)))*A2
  add(oo,0,('odd',p,c))
 from math import factorial,lcm
 for i in range(15):
  j=14-i; terms={}
  norm=F(1)
  for r in range(8):
   if r:norm*=-F(100,2*r*(96+28-2*r-2))
   for u in range(r+1):
    v=r-u
    if i>=2*u and j>=2*v:terms[i-2*u,j-2*v]=norm*F(factorial(r)*factorial(i)*factorial(j),factorial(u)*factorial(v)*factorial(i-2*u)*factorial(j-2*v))
  L=lcm(*(x.denominator for x in terms.values()))
  add(lambda p,c,s:sum(int(L*a)*p**u*c**v for (u,v),a in terms.items()),0,('harm14',i,j))
 return rows,vals,names,B,forbid

if __name__=='__main__':
 import json
 from pathlib import Path
 R,V,M,B,Fb=matrix()
 assert len(B)==222 and not Fb and len(R)==93
 D={int(k):F(v) for k,v in json.loads(Path(__file__).with_name('anchor-counts.json').read_text()).items()}
 assert len(D)==45 and all(v>0 for v in D.values())
 assert D[0]==1 and B[0]==((-10,0,-10),(0,-10,10),(0,10,10),(10,0,-10))
 for r,w,m in zip(R,V,M):
  assert sum(r[i]*v for i,v in D.items())==w,m
 print('PASS 222 potential orbits, 45 rational support orbits, 93 marked constraints.')
