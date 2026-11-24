from fractions import Fraction as Q
from math import comb
from functools import lru_cache
N=22
def mul(a,b):
 return tuple(sum(a[j]*b[i-j] for j in range(i+1)) for i in range(N+1))
def scale(a,s):return tuple(a0*s for a0 in a)
def add(a,b): return tuple(a0+b0 for a0,b0 in zip(a,b))
def sub(a,b):return add(a,scale(b,-1))
@lru_cache(None)
def pw(a,n):
 if not n:return (1,)+(0,)*N
 return mul(a,pw(a,n-1))
def eta(data):
 num=sum(d*k for d,k in data.items())
 assert num%24==0 and num>=0
 sh=num//24; a=(0,)*sh+(1,)+(0,)*(N-sh)
 for d,k in data.items():
  for j in range(d,N+1,d):
   b=[0]*(N+1)
   for i in range(min(N//j,k)+1): b[i*j]=(-1)**i*comb(k,i)
   a=mul(a,b)
 return a
def sigma(n,k):return sum(i**k for i in range(1,n+1) if n%i==0)
def E(d,k):
 return (1,)+tuple(({2:-24,4:240}[k]*sigma(i//d,k-1) if i%d==0 else 0) for i in range(1,N+1))
Ds=(1,2,5,10)
e=tuple(Q(sum(d*d*E(d,4)[j] for d in Ds),130) for j in range(N+1))
phi=tuple(Q(sum(sign*d*E(d,2)[j] for d,sign in zip(Ds,(-1,-1,1,1))),12) for j in range(N+1))
f=scale(sub(pw(phi,2),e),Q(13,28))
assert f==add(eta({1:4,5:4}),scale(eta({2:4,10:4}),4))
g=eta({1:4,2:4,5:4,10:4})
assert e[0]==1 and f[:2]==(0,1) and g[:4]==(0,0,0,1)
expected={48:157248000,64:-79877913600,72:-125084010240,80:-81289105920000,88:-276587713406400,96:-48713592107980800}
for n,want in expected.items():
 h=n//8;m=3*n//16;B={}
 for b in range(h//2+1):
  for a in range(h-2*b+1):
   c=h-a-2*b
   B[a+3*b]=mul(mul(pw(e,c),pw(f,a)),pw(g,b))
 assert set(B)==set(range(m+1))
 v=B[0]
 for j in range(1,m+1):v=sub(v,scale(B[j],v[j]))
 assert v[:m+1]==(1,)+(0,)*m and v[m+2]==want
 if n==48:assert v[m+1]==144668160
 if want<0:
  # The norm-2m minimum branch differs only by A times the last cusp basis.
  coeff=B[m][m+2]
  if coeff>0: print(n,"squared min <=",2*m,"necessary norm",2*m,"count >=",Q(-want,coeff))
  else: print(n,"extremal norm",2*m+2,"excluded")
for n,m,h,want in [(16,4,4,8),(32,6,6,18),(48,8,4,Q(162,5)),(64,10,6,50),(80,12,8,72),(96,14,6,Q(512,5))]:
 var=Q(5*m*m,n//2+1)
 assert var > (16 if n==80 else 4 if n in (32,64,96) else 0)
 assert Q((m-2)**2*5,8)<want
 assert max(Q(m*m,2),Q((5*m-h)**2,40))==want
 assert m>Q(n//2+6,4)
assert Q(162,5)<36
print("PASS level-10 Fourier obstructions and full golden twist bounds on even 2-modular parents")
