"""Exact leading Fourier coefficient ranks for selected golden Hilbert weights."""
from fractions import Fraction as F
from itertools import permutations,product
from collections import defaultdict
from math import isqrt,gcd,comb
# Use rational pairs x+y sqrt(5).
Z=(F(0),F(0)); O=(F(1),F(0))
def add(x,y): return x[0]+y[0],x[1]+y[1]
def neg(x): return -x[0],-x[1]
def mulq(x,y): return x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def div(x,y):
 d=y[0]**2-5*y[1]**2
 return mulq(x,(y[0]/d,-y[1]/d))
def pw(x,n):
 z=O
 for _ in range(n): z=mulq(z,x)
 return z
T=8
def inds(T):
 return [(0,0)]+[(n,b) for n in range(1,T+1) for b in range(-isqrt(5*n*n),isqrt(5*n*n)+1) if (n-b)%2==0]
def factors(x):
 for p in range(2,isqrt(x)+1):
  if p*p>x:break
  k=0
  while x%p==0:k+=1;x//=p
  if k:yield p,k
 if x>1:yield x,1
def idsum(n,b,k):
 N=(5*n*n-b*b)//4
 G=gcd(n,(b-n)//2)
 A=1
 for p,e in factors(N):
  if p==5:L=[(5,e)]
  elif p%5 in (2,3):
   assert e%2==0
   L=[(p*p,e//2)]
  else:
   t=0
   while G%p==0:G//=p;t+=1
   assert t<=e-t
   L=[(p,t),(p,e-t)]
  for a,b in L:A*=sum(a**((k-1)*i) for i in range(b+1))
 return A
def eis(k,c):
 return {(0,0):O}|{(n,b):(c*idsum(n,b,k),F(0)) for n,b in inds(T)[1:]}
def times(x,y):
 z=defaultdict(lambda:Z)
 for (n,b),a in x.items():
  for (m,c),d in y.items():
   if n+m<=T:
    j=(n+m,b+c)
    z[j]=add(z[j],mulq(a,d))
 return {j:v for j,v in z.items() if v!=Z}
def linear(*fs):
 z=defaultdict(lambda:Z)
 for a,x in fs:
  for j,v in x.items():z[j]=add(z[j],mulq(a,v))
 return {j:v for j,v in z.items() if v!=Z}
def powf(x,n):
 z={(0,0):O}
 for _ in range(n):z=times(z,x)
 return z
def normalize(f,j):
 c=div(O,f[j])
 return linear((c,f))
BN=[F(1)]
for k in range(1,11):
 BN.append(-sum(comb(k+1,j)*BN[j] for j in range(k))/(k+1))
def Bern(k,x):return sum(comb(k,j)*BN[j]*x**(k-j) for j in range(k+1))
def escale(k):
 zeta=-BN[k]/k
 L= -F(5**(k-1),k)*sum(ch*Bern(k,F(j,5)) for j,ch in enumerate((1,-1,-1,1),1))
 return F(4)/(zeta*L)
assert [escale(k) for k in (2,6,10)]==[F(120),F(2520,67),F(6600,412751)]
E,E6,E10=(eis(k,escale(k)) for k in (2,6,10))
C6=normalize(linear((O,powf(E,3)),(neg(O),E6)),(1,-1))
C10=linear((O,powf(E,5)),(neg(O),E10))
E2C6=times(powf(E,2),C6)
C10=normalize(linear((O,C10),(neg(div(C10[(1,-1)],E2C6[(1,-1)])),E2C6)),(2,-2))
# The square root normalized to coefficient 1 at (1,-1).
S5={}
for n,b in inds(T-1)[1:]:
 j=(n+1,b-1)
 a=O if j==(2,-2) else div(add(C10.get(j,Z),neg(times(S5,S5).get(j,Z))),(F(2),F(0)))
 if a!=Z:S5[n,b]=a
 assert times(S5,S5).get(j,Z)==C10.get(j,Z)
def der(f,e,side):
 return {(n,b):mulq(a,pw((F(5*n),F(side*b)),e)) for (n,b),a in f.items()}
def rc(f,g,k,l,m):
 return linear(*[((F((-1)**(m-a)*comb(m+l-1,a)*comb(m+k-1,m-a)),F(0)),
                 times(der(f,a,-1),der(g,m-a,-1))) for a in range(m+1)])
# Gundlach's additional weight-15 generator via the modular Jacobian, divided by S5.
fs=[E,C6,C10]; ws=[2,6,10]
def D(f,axis):return {(n,b):mulq(a,(F(n if axis==0 else b),F(0))) for (n,b),a in f.items() if n}
J=linear(*[((F(ws[i]*(-1)**((i>j)+(i>k)+(j>k))),F(0)),
           times(times(fs[i],D(fs[j],0)),D(fs[k],1))) for i,j,k in permutations(range(3))])
S15={}
for n,b in inds(T-1)[1:]:
 j=(n+1,b-1)
 a=add(J.get(j,Z),neg(times(S15,S5).get(j,Z)))
 if a!=Z:S15[n,b]=a
 assert times(S15,S5).get(j,Z)==J.get(j,Z)
assert S15 and S5
# R = C[E,S5,C6], symmetric for the signed exchange convention.
def rbasis(k):
 for a,b,c in product(range(k//2+1),range(k//5+1),range(k//6+1)):
  if 2*a+5*b+6*c==k:yield times(times(powf(E,a),powf(S5,b)),powf(C6,c))
G={
 0:[(0,{(0,0):O}),(15,S15)],
 2:[(7,rc(E,S5,2,5,1)),(8,rc(E,C6,2,6,1)),(11,rc(S5,C6,5,6,1))],
 4:[(4,rc(E,E,2,2,2)),(7,rc(E,S5,2,5,2))]
}
def rank(M):
 M=[list(r) for r in M]
 R=0
 for c in range(len(M[0])):
  j=next((j for j in range(R,len(M)) if M[j][c]!=Z),None)
  if j is None:continue
  M[R],M[j]=M[j],M[R]
  for i in range(R+1,len(M)):
   if M[i][c]==Z:continue
   u=div(M[i][c],M[R][c])
   for h in range(c,len(M[0])):M[i][h]=add(M[i][h],neg(mulq(u,M[R][h])))
  R+=1
 return R
bad=[j for j in inds(4) if j[0]<4 or j[1]!=0]

# Scalar Hilbert weight 24, both signs, including the antisymmetric generator.
from replay import size
C=[times(g,f) for w,g in G[0] for f in rbasis(24-w)]
rows=[[(g.get(t,Z)) for g in C] for t in bad+[(4,0)]]
# Solve simultaneously the constant and norm-8 coefficient parts.
M=[list(r)+[(F(t==(0,0)),F(0)),(F(t==(4,0)),F(0))]
   for r,t in zip(rows,bad+[(4,0)])]
N=len(C)
for i in range(N):
    j=next(j for j in range(i,len(M)) if M[j][i]!=Z)
    M[i],M[j]=M[j],M[i]
    M[i]=[div(v,M[i][i]) for v in M[i]]
    for j in range(len(M)):
        if j==i:continue
        c=M[j][i]
        M[j]=[add(v,neg(mulq(c,w))) for v,w in zip(M[j],M[i])]
assert all(v==Z for row in M[N:] for v in row)
for n,b in inds(8):
    if n<4:continue
    for a,A in enumerate((0,1)):
        x=Z
        for g,row in zip(C,M):
            x=add(x,mulq(g.get((n,b),Z),row[N+a]))
        assert x==(F(size(A,2*n,abs(b))-(size(0,2*n,abs(b)) if a else 0)),F(0)),(n,b,a,x)
print('PASS exact scalar Hilbert theta coefficients through squared norm 16')
