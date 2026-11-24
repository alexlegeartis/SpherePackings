"""Exact two-shell propagation certificate for the anisotropic double-radical class."""
from itertools import combinations,product
from collections import Counter
from fractions import Fraction as F
from math import isqrt
from subprocess import run

G=[[0]*6 for _ in range(6)]
for i in range(6):G[i][i]=2
for i,j in [(0,1),(1,2),(2,3),(3,4),(2,5)]:G[i][j]=G[j][i]=-1
R=(2,1,0,2,1,0)
half_inv_diag=(4,10,18,10,4,6)  # diagonal of 3*G^{-1}
assert all(sum(G[i][j]*R[j] for j in range(6))%3==0 for i in range(6))
assert {v for v in product(range(3),repeat=6) if all(sum(G[i][j]*v[j] for j in range(6))%3==0 for i in range(6))} == {tuple(k*x%3 for x in R) for k in range(3)}
def Q(p):return sum(p[i]*p[j]*G[i][j] for i in range(6) for j in range(6))
assert all(G[i][i]%2==0 for i in range(6))
def normtype(w):
 return Q(w)//2%3
def canon(w):return min(tuple((a+k*b)%3 for a,b in zip(w,R)) for k in range(3))
def sd(w,v):return tuple((a-b)%3 for a,b in zip(w,v))
L=sorted(set(map(canon,product(range(3),repeat=6))))
I={v:i for i,v in enumerate(L)}
assert len(L)==243
assert [sum(normtype(x)==i for x in L) for i in range(3)]==[81,72,90]
# Actual short lifts suffice, so no assertion about global minima is needed.
# E6 roots' reductions are all of type one and their quotient classes are distinct.
# The shell enumeration below directly verifies each relevant type.
B=[isqrt(6*a//3) for a in half_inv_diag]
depth={}
for p in product(*(range(-a,a+1) for a in B)):
 n=Q(p)
 if n<=6:
  v=canon(p);depth[v]=min(depth.get(v,n),n)
B8=[isqrt(8*a//3) for a in half_inv_diag]
short8=set()
for p in product(*(range(-a,a+1) for a in B8)):
 if Q(p)<=8:short8.add(tuple(x%3 for x in p))
assert len(short8)==727
assert set(product(range(3),repeat=6))-short8 == {tuple(k*x%3 for x in R) for k in (1,2)}
assert len(depth)==243
assert all(depth[v]==2*normtype(v) if normtype(v)!=0 else depth[v]==(0 if v==L[0] else 6) for v in L)
T=[[normtype(sd(v,w)) for w in L] for v in L]
nb=[sum(1<<j for j,z in enumerate(row) if z==1) for row in T]
N=243

def cliques(P,mask,k):
 if not k: yield P; return
 while mask.bit_count()>=k:
  b=mask&-mask;mask^=b;i=b.bit_length()-1
  yield from cliques(P+(i,),mask&nb[i],k-1)

roots=set(cliques((0,),nb[0],5))
assert len(roots)==432
S=min(roots)
maps=[]
for j in range(6):
 ref=[]
 for v in L:
  w=list(v);w[j]-=sum(G[j][i]*v[i] for i in range(6));ref.append(I[canon(w)])
 maps.append(ref)
maps.append([I[canon(tuple(-a for a in v))] for v in L])
seen={S};todo=[S]
for P in todo:
 for U in maps:
  Z=tuple(sorted(U[i] for i in P))
  if Z not in seen:seen.add(Z);todo.append(Z)
assert seen==roots
assert not list(cliques((0,),nb[0],6))
assert [sum(sum(T[v][s]==1 for s in S)==k for v in range(N)) for k in range(7)]==[0,132,75,0,30,6,0]


# Every affine K6 is equivalent to S (translate one vertex to zero).
W0=sorted(set(S)|{v for v in range(N) if sum(T[s][v]==1 for s in S)>=4})
assert len(W0)==36
D0=W0+sorted(set(range(N))-set(W0))
pos={v:i for i,v in enumerate(W0)}
bits=lambda L:sum(1<<pos[v] for v in L)
mat=[]
for v in D0:
 R0=[u for u in W0 if T[v][u]==1]
 masks=[bits(P) for P in combinations(R0,4) if all(T[x][y]!=0 for x,y in combinations(P,2))]
 mat.extend([len(masks),*masks])
def certify(om,h):
 seed=set(S)-set(om)
 raw=[N,h,len(W0),150,int(not om),bits(seed),bits(set(W0)-set(S)),bits(om)]
 run(['./closure'],input=' '.join(map(str,raw+mat)),text=True,check=True)
assert F(48*(5*6+30),243)<12
certify((),11)
for t in S:certify((t,),6)
for t in combinations(S,2):certify(t,1)
assert 150-48 > 9*9
target=F(64,9)*(F(9,5)*F(18,5))**3
assert 3*(243+2*194)==1893<target
print('PASS every affine vertical-nontrivial radical invariant motif falls short')
