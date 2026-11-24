"""First Fourier-Jacobi coefficients of even unimodular rank-96 theta."""
from fractions import Fraction as F
from collections import defaultdict
from itertools import product
from functools import lru_cache
import sympy as sy
T=10
ZERO={};ONE={(0,0):1}
def mul(a,b):
 d=defaultdict(int)
 for (n,r),u in a.items():
  for (m,s),v in b.items():
   if n+m<=T:d[(n+m,r+s)]+=u*v
 return {k:v for k,v in d.items() if v}
def powd(a,p):
 b=ONE
 for _ in range(p):b=mul(a,b)
 return b
def add(*terms):
 v=defaultdict(int)
 for c,a in terms:
  for j,k in a.items():v[j]+=c*k
 return {j:k for j,k in v.items() if k}
def eis(w,c):
 return {(0,0):1}|{(n,0):c*sum(d**(w-1) for d in range(1,n+1) if n%d==0) for n in range(1,T+1)}
E2,E4,E6=eis(2,-24),eis(4,240),eis(6,-504)
D={j:v//1728 for j,v in add((1,powd(E4,3)),(-1,powd(E6,2))).items()}
assert D[(1,0)]==1
J={(0,-1):1,(0,0):-2,(0,1):1}
from math import comb
for n in range(1,T+1):
 J=mul(J,mul(mul(powd({(0,0):1,(n,1):-1},2),powd({(0,0):1,(n,-1):-1},2)),{(n*k,0):comb(k+3,3) for k in range(T//n+1)}))
K=add((1,{(n,r):(-24*n+6*r*r)*v for (n,r),v in J.items()}),(-5,mul(E2,J)))
assert all(K.get((0,r))==v for r,v in ((-1,1),(0,10),(1,1)))
@lru_cache(None)
def polys(w):
 return [mul(powd(D,k),mul(powd(E4,a),powd(E6,b))) for k in range(w//12+1) for a in range(w//4+1) for b in range(2) if 12*k+4*a+6*b==w]
def ebasis(m,w=48):
 return [mul(mul(powd(J,a),powd(K,m-a)),E) for a in range(m+1) for E in polys(w+2*a)]
def theta(A,w=48):
 E=polys(w)
 def coeff(a,n):return a.get((n,0),0)
 x=list(sy.Matrix([[coeff(a,n) for a in E] for n in range(5)]).inv()*sy.Matrix([1,0,0,0,A]))
 return {n:sum(t*coeff(a,n) for t,a in zip(x,E)) for n in range(T+1)}
def rank2neg(m,n,r,mi):
 D=4*m*n-r*r
 if D<0:return True
 if D==0:return False
 # Binary quadratic form via Gauss reduction
 x,y,z=m,n,r
 while True:
  if x>y:x,y=y,x
  if abs(z)<=x:break
  l=(z+x)//(2*x)
  y,z=y-l*z+l*l*x,z-2*l*x
 return min(x,y)<mi
def row(E,n,r):return [a.get((n,r),0) for a in E]
def mj(m,mi,ta,other={}):
 E=ebasis(m);N=ta[m];Q=[];Z=[]
 def use(L,b):Q.append(L);Z.append(b)
 def get(n,r):return row(E,n,r)
 for n in range(T+1):
  # z=0 sums
  use([sum(v for (i,r),v in a.items() if i==n) for a in E],N*ta[n])
  # no short vector pairs outside admissible support
  for r in range(-m,m+1):
   if r==0 and n==0:use(get(n,r),N);continue
   if rank2neg(m,n,r,mi) or n<mi and r:use(get(n,r),0)
 for n,otherE in other.items():
  for r,val in otherE.items():use(get(n,r),val)
 mat=sy.polys.matrices.DomainMatrix.from_Matrix(sy.Matrix([l+[b] for l,b in zip(Q,Z)]))
 rr,piv=mat.rref()
 R=rr.to_Matrix();s=len(E)
 assert s not in piv, (m,mi,piv)
 Fv=sorted(set(range(s))-set(piv))
 Bs=[]
 for g in [None]+Fv:
  x=[int(j==g) for j in range(s)]
  for i,j in enumerate(piv): x[j]=R[i,s] if g is None else -R[i,g]
  Bs.append(x)
 def coeff(n,r):return tuple(sum(u*v for u,v in zip(get(n,r),x)) for x in Bs)
 return coeff, len(Fv)
if __name__=='__main__':
 import json
 from pathlib import Path
 t=theta(0)
 j5,d5=mj(5,5,t)
 assert d5==0 and len(ebasis(5))==29
 j6,d6=mj(6,5,t,{5:{r:j5(6,r)[0] for r in range(-12,13)}})
 assert d6==0 and len(ebasis(6))==35
 assert j5(6,5)[0]==j6(6,7)[0]==11403180502176153600000
 C=json.loads(Path(__file__).with_name('round17-counts.json').read_text())
 N,W=C['min']
 prev={p:F(N)*sum(F(w) for r,b,w in W if r==p)/(2 if p else 1) for p in range(8)}
 assert all(prev[p]<=j6(6,p)[0] for p in range(7))
 assert prev[7]>F(1538,1000)*j6(6,7)[0]
 print('PASS extremal norm-12 pair budget at fixed dot +7:',j6(6,7)[0],'; old witness/whole shell =',float(prev[7]/j6(6,7)[0]))
