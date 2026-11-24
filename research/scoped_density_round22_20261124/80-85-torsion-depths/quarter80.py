from fractions import Fraction as F
from math import ceil,floor,isqrt
from sympy import Matrix
N=16
def C(x):return (x,)+(0,)*(N-1)
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,v):return tuple(v*x for x in a)
def mul(a,b):return tuple(sum(a[i]*b[j-i]for i in range(j+1))for j in range(N))
def p(a,n):
 s=C(1)
 for _ in range(n):s=mul(s,a)
 return s
th3=tuple(1 if k==0 else (2 if isqrt(k)**2==k else 0) for k in range(N))
th4=tuple(1 if k==0 else (2*(-1)**isqrt(k) if isqrt(k)**2==k else 0) for k in range(N))
th2=tuple(2*any(k==n*(n+1) for n in range(N))for k in range(N))
X=(0,)+p(th2,4)[:-1];Y=p(th4,4);W=add(X,Y)
def quarter(r,H,E):
 alpha=F(4-r,4) if r else F(0);J=ceil(H/2-alpha)
 i=ceil(2*(E-F(r,8)));d=20-bool(r)-J-i
 if d<0:return None
 def seq(X,Y,W,P):
  return [mul(P,mul(p(Y,i),mul(p(W,J),mul(p(X,u),p(Y,d-u)))))for u in range(d+1)]
 P=mul(p(th4,r),p(th3,4-r)) if r else C(1)
 A=seq(X,Y,W,P)
 try:a,v=Matrix([[g[t]for g in A]for t in range(4)]).gauss_jordan_solve(Matrix([1,0,0,0]))
 except ValueError:return None
 J0=seq(W,scale(Y,-1),X,mul(p(th4,r),p(th2,4-r)) if r else C(1))
 V=tuple(sum(x*g[t]for x,g in zip(a,J0))for t in range(8))
 return i,J,d,2*(alpha+J),V,tuple(v),2**(r+4*i-40)*a[-1]
def geom(r,H):
 parity=F((r-int(2*H))%4,4);assert parity in (0,F(1,2))
 m=parity
 while m*m<F(7,20)*H:m+=1
 q=max(F(13,2)+H/4-parity/2,7+H/4-m/2)
 return F(r,8)+F(floor(2*(q-F(r,8))),2)
best=F(7);kept=[]
for r in range(4):
 for kk in range(4,15):
  H=F(kk,2)
  if (H-F(r,2)).denominator>1:continue
  G=geom(r,H)
  if G<=7:continue
  D=F(r,8)+F(floor(2*(7-F(r,8)))+1,2)
  while D<=G:
   B=quarter(r,H,D)
   if B is not None:
    print('candidate',r,H,'geom',G,'depth at least',D,'data',B,flush=True)
    best=max(best,D);kept.append((r,H,D))
   D+=F(1,2)
print('best mod gate',best,kept)

G=[(r,H,d)for r,H,d in kept if d>F(29,4)]
assert G==[(0,F(4),F(15,2)),(3,F(7,2),F(59,8)),(3,F(9,2),F(59,8))]
for r,H,d in G:
 i,j,k,s,v,free,c=quarter(r,H,d)
 assert not free and 0<s
 want=abs(v[j])
 bound=160 if r==0 else 2*80*(1-F(1,81))/(1-F(80,81))
 assert (r,s,want)==((0,F(4),9920)if r==0 else (3,F(9,2),17572))
 assert want>bound
 print('SIGNED SHELL EXCLUDED',r,H,d,'required',want,'max',bound)
print('PASS universal fourth-torsion squared depth <= 29/4 for rank80 min-eight L')
