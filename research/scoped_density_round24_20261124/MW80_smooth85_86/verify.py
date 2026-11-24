from collections import Counter
from itertools import product
import json
p=41
O=None
def add(x,y):
 return ((x%p+y%p)%p)+p*((x//p+y//p)%p)
def neg(x):
 return (-x%p)+p*(-(x//p)%p)
def sub(x,y):return add(x,neg(y))
def mul(x,y):
 a,b=divmod(x,p);c,d=divmod(y,p)
 return (b*d+3*a*c)%p+p*((a*d+b*c)%p)
def pw(x,n):
 if n<0:n%=p*p-1
 y=1
 while n:
  if n&1:y=mul(y,x)
  x=mul(x,x);n>>=1
 return y
def group(A,B):
 if A is O:return B
 if B is O:return A
 x,y=A;X,Y=B
 if x==X and add(y,Y)==0:return O
 if x==X:m=mul(mul(3,mul(x,x)),pw(mul(2,y),-1))
 else:m=mul(sub(Y,y),pw(sub(X,x),-1))
 z=sub(sub(mul(m,m),x),X)
 return z,sub(mul(m,sub(x,z)),y)
def times(P,n):
 if n<0:
  return times((P[0],neg(P[1])), -n)
 Q=O
 while n:
  if n&1:Q=group(P,Q)
  P=group(P,P)
  n>>=1
 return Q
def polmul(X,Y):
 w=[0]*(len(X)+len(Y)-1)
 for j,x in enumerate(X):
  for k,y in enumerate(Y):w[j+k]=add(w[j+k],mul(x,y))
 return w
def polpow(P,n):
 out=[1]
 for j in range(n):out=polmul(out,P)
 return out
from fractions import Fraction as F
from math import factorial
P0,P1=(3,328),(50,14)
def model(P):
    x,y=P
    return add(pw(x,3),1)==pw(y,2)
assert pow(3,20,41)==40
assert model(P0) and model(P1)
assert times(P0,42) is None and times(P1,42) is None
prim0=[times(P0,q)for q in range(42)]
prim1=[times(P1,q)for q in range(42)]
labels={group(a,b):(i,j)for i,a in enumerate(prim0) for j,b in enumerate(prim1)}
assert len(labels)==1764
C=Counter(mul(y,y)for y in range(p*p))
assert len(labels)==1+sum(C[add(pw(x,3),1)]for x in range(p*p))

# All elementary polynomial norm-14 sections are compressed by their
# constant terms. Any last column (b,d) of norm 1 can be completed
# with (a,c)=(-d^p,b^p), preserving t^42+1.
X={neg(pw(b,14)):b for b in range(p*p)}
Y={pw(d,21):d for d in range(p*p)}
stock={}
for x,y in product(X,Y):
    if not model((x,y)): continue
    b,d=X[x],Y[y]
    a,c=neg(pw(d,41)),pw(b,41)
    assert add(pw(a,42),pw(c,42))==1 and add(pw(b,42),pw(d,42))==1
    assert add(mul(pw(a,41),b),mul(pw(c,41),d))==0
    stock[labels[(x,y)]]=(a,b,c,d)
assert len(stock)==239
for q in [2,3,7]:
    s=42//q
    z=[(a//s,b//s)for a,b in stock if a%s==b%s==0]
    assert any((a*d-b*c)%q for (a,b),(c,d) in product(z,z))
for t in [1,2]:
    assert any(l[t-1] for l in stock)

def color(x,M):
    return tuple(sum(y*row[k]for y,row in zip(x,M))%42 for k in (0,1))
D=[(25,23),(0,13),(21,37),(2,16),(9,13)]
E=[(8,18),(12,1),(15,1),(25,21),(0,2),(9,23)]
adj=[(0,1),(1,2),(2,3),(3,4),(2,5)]
contacts={}
counts={}
cnt=0
for x in product(range(-3,4),repeat=5):
    n=sum(y*y for y in x)
    if n==0 or n>=16 or sum(x)%2:continue
    a,b,c=tuple(sum(x[:k])for k in range(1,4))
    d,e=(c+x[3]-x[4])//2,(c+x[3]+x[4])//2
    z=(a,b,c,d,e); assert (a,-a+b,-b+c,-c+d+e,-d+e)==x
    xy=color(z,D); assert any(xy);cnt+=1
    if n==2 and xy in stock:contacts[85]=(z,xy)
assert cnt==2330
counts[85]=cnt

# E6 Cartan: positive leading principal minors and exact inverse diagonal.
def det(A):
    B=[[F(x)for x in r]for r in A];n=len(B);p=F(1)
    for j in range(n):
        hj=next((k for k in range(j,n)if B[k][j]),None)
        if hj is None:return F(0)
        if hj!=j:B[hj],B[j]=B[j],B[hj];p=-p
        z=B[j][j];p*=z
        for i in range(j+1,n):
            y=B[i][j]/z
            for k in range(j+1,n):B[i][k]-=y*B[j][k]
    return p
S=[[2 if i==j else (-1 if (min(i,j),max(i,j)) in adj else 0)for j in range(6)]for i in range(6)]
assert all(det([r[:i]for r in S[:i]])>0 for i in range(1,7))
assert det(S)==3
from math import isqrt
diag=[det([[v for j,v in enumerate(row)if j!=i]for k,row in enumerate(S)if k!=i])/3 for i in range(6)]
W=[isqrt((14*d).__floor__())for d in diag]
assert W==[4,6,9,6,4,5]
cnt=0
for z in product(*[range(-r,r+1)for r in W]):
    n=2*(sum(y*y for y in z)-sum(z[i]*z[j]for i,j in adj))
    if n==0 or n>=16:continue
    assert n>=2
    xy=color(z,E);assert any(xy);cnt+=1
    if n==2 and xy in stock:contacts[86]=(z,xy)
assert cnt==9972
counts[86]=cnt

pi=16*sum(F((-1)**j,(2*j+1)*5**(2*j+1))for j in range(40))-4*sum(F((-1)**j,(2*j+1)*239**(2*j+1))for j in range(41))
for d,a,c in [(85,F('3.3809404126266e-19'),4),(86,F('2.1043362358080e-19'),3)]:
    v2=(2**86*pi**42*F(factorial(43),factorial(86)))**2 if d==85 else (pi**43/factorial(43))**2
    assert v2*F(4**d,c*(41**6*1764)**2)>a*a
    assert d in contacts
    z,v=contacts[d];a0,b0,c0,d0=stock[v]
    x=[neg(u)for u in polpow([b0,a0],14)]
    y=polpow([d0,c0],21)
    f=polpow(x,3);f[0]=add(f[0],1);f[-1]=add(f[-1],1)
    assert polpow(y,2)==f
    print(d,'squared min 16',counts[d],'full vertical kernel exclusions','contact',z,v,stock[v])
