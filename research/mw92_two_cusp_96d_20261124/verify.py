"""Replay the two-cusp test in F_47[i], with no numerical search."""
from fractions import Fraction as F
from collections import Counter
from pathlib import Path
import json
p=47; O=(0,0); I=(1,0)
def a(x,y):return tuple((u+v)%p for u,v in zip(x,y))
def m(x,y):
 u,v=x;r,s=y
 return ((u*r-v*s)%p,(u*s+v*r)%p)
def neg(x):return tuple(-z%p for z in x)
def P(x,k):
 y=I
 while k:
  if k&1:y=m(y,x)
  x=m(x,x);k//=2
 return y
def inv(x):assert x!=O;return P(x,p*p-2)
def canon(A):
 A=list(A)
 while A and A[-1]==O:A.pop()
 return A
def add(A,B):
 C=[O]*max(len(A),len(B))
 for j,x in enumerate(A):C[j]=a(C[j],x)
 for j,x in enumerate(B):C[j]=a(C[j],x)
 return canon(C)
def sc(A,c):return canon([m(z,c)for z in A])
def sub(A,B):return add(A,sc(B,(-1,0)))
def mul(A,B):
 if not A or not B:return []
 C=[O]*(len(A)+len(B)-1)
 for i,u in enumerate(A):
  for j,v in enumerate(B):C[i+j]=a(C[i+j],m(u,v))
 return canon(C)
def power(A,k):
 B=[I]
 for j in range(k):B=mul(B,A)
 return B
def quo(A,B):
 assert B
 A=canon(A);C=[O]*max(0,len(A)-len(B)+1);z=inv(B[-1])
 while len(A)>=len(B):
  u=m(A[-1],z);k=len(A)-len(B);C[k]=u
  A=sub(A,[O]*k+sc(B,u))
 return canon(C),A
def gcd(A,B):
 while B:A,B=B,quo(A,B)[1]
 return sc(A,inv(A[-1]))
def decode(A):return [(v%p,v//p)for v in A]
def evaluate(A,t):
 y=O
 for z in reversed(A):y=a(m(y,t),z)
 return y
HERE=Path(__file__).resolve().parent
f=HERE.parent/'mw92_derivative_lines_96d_20261124'/'witness_polys.json'
if not f.exists():f=HERE/'witness_polys.json'
J=json.loads(f.read_text())
def get(n):
 A,B=[decode(J[n][x])for x in ('X','Y')]
 assert (len(A),len(B))==(19,28)
 assert sub(power(B,2),power(A,3))==[O,(-1%p,0)]+[O]*45+[I]
 return A,B
A,B=get('2');h=(21,17);assert P(h,6)==(21,0)
L=[(23,0),I];N=[(38,0),(23,0)]
def trans(A,d):
 H=[]
 for j,z in enumerate(A):H=add(H,sc(mul(power(N,j),power(L,d-j)),z))
 return H
C=sc(trans(A,18),P(inv(h),2));D=sc(trans(B,27),P(inv(h),3))
X=mul(A,power(L,2));Y=mul(B,power(L,3))
J0=sub(X,C);J1=sub(Y,D)
assert gcd(J0,J1)==[(20,0),(25,0),I]
fnum=sub(power(add(Y,D),2),mul(add(X,C),power(J0,2)))
fden=mul(power(L,2),power(J0,2))
G=gcd(fnum,fden);fnum,e=quo(fnum,G);fden,e1=quo(fden,G)
assert not e and not e1 and len(fnum)==21
assert sc(fden,inv(fden[-1]))==[(24,0),(13,0),(7,0),(3,0),I]
assert mul([(-28%p,0),I],[(-41%p,0),I])==[(20,0),(25,0),I]
print('PASS: explicit difference with reduced poles at exactly 28 and 41')
A,B=get('0');g=(2,1);N=p*p-1
G=[P(g,i)for i in range(N)];assert len(set(G))==N
log={u:j for j,u in enumerate(G)}
w=[m(evaluate(A,(b,0)),inv(evaluate(B,(b,0)))) for b in range(p)]
assert O not in w
C=Counter(log[z]%8 for z in w)
assert [C[j]for j in range(8)]==[7,6,3,7,8,4,6,6]
H=[z for z in G if P(z,6)[1]==0];assert set(H)==set(G[::8])
assert len({m(r,z)for r in w for z in H})==N
print('PASS: every nonzero second-cusp colour occurs at norm 18')
t=F(18*p,p-1);pi=F(314,100)
assert 64*F(t/18)**48/(pi**2*t**2)<F(116,1000)
assert F(16,p*p)*F(10,9)**48/pi**2<F(116,1000)
print('PASS: density ratio at both critical endpoints < 11.6%')
