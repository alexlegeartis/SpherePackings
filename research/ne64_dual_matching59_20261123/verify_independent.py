#!/usr/bin/env python3
"""Independent integer-only arithmetic replay of the literal period and motif."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import json
assert __debug__
P=Path(__file__).resolve().parent
def read(p):return [list(map(int,l.split()))for l in (P/p).read_text().splitlines()]
def tr(x):return [list(r)for r in zip(*x)]
def m(a,b):return [[sum(x*y for x,y in zip(v,w))for w in tr(b)]for v in a]
def eye(n):return [[int(i==j)for j in range(n)]for i in range(n)]
def scale(a,t):return [[t*x for x in v]for v in a]
def sub(a,b):return [[x-y for x,y in zip(v,w)]for v,w in zip(a,b)]
def sel(a,I):return [[v[i]for i in I]for v in a]
def det(G, positive=False):
 a=[v[:]for v in G];n=len(a);p=1;sgn=1
 for k in range(n):
  j=next((j for j in range(k,n)if a[j][k]),None)
  if j is None:return 0
  if j!=k:a[j],a[k]=a[k],a[j];sgn=-sgn
  z=a[k][k]
  if positive:assert j==k and z>0
  for i in range(k+1,n):
   for j in range(k+1,n):
    x,y=divmod(z*a[i][j]-a[i][k]*a[k][j],p);assert not y;a[i][j]=x
   a[i][k]=0
  p=z
 return sgn*p
T=read('change.txt');assert abs(det(T))==1
rows=read('ne64.triangle');assert list(map(len,rows))==list(range(1,65))
G=[[rows[max(i,j)][min(i,j)]for j in range(64)]for i in range(64)]
assert det(G,True)==3**32
C=m(m(T,G),tr(T));J=read('data/ambient_gram.txt');assert J==tr(J)
assert m(C,J)==scale(eye(64),3)
A=read('normal.txt')
B=[[12 if i==j else (-4 if {i,j} in ({1,2},{3,4}) else -2)for j in range(5)]for i in range(5)]
assert m(m(A,C),tr(A))==B
I=[0,1,3,5,25];free=sorted(set(range(64))-set(I))
assert abs(det(sel(A,I)))==1
N=read('data/period_basis.txt');Q=read('data/period_gram.txt')
assert sel(N,free)==eye(59) and m(N,tr(A))==[[0]*5 for _ in range(59)]
assert m(m(N,J),tr(N))==Q and det(Q,True)==32768*3**28
L=read('data/lifts.txt');S=read('data/shift_numerators.txt')
assert m(L,tr(A))==[[0]*5]+eye(5)
Y=m(S,N);assert m(Y,tr(A))==[[0]*5 for _ in range(6)]
M=sub(scale(L,96),Y)
Jm=m(M,J)
# Find integral coefficients for the normal projection. The five pivots
# are unimodular, so each row of Jm has a unique coefficient vector.
Ai=sel(A,I);K=det(Ai)
Vi=[[((-1)**(i+j))*det(sel([r for z,r in enumerate(Ai)if z!=j],[u for u in range(5)if u!=i]))//K
     for j in range(5)]for i in range(5)]
assert m(Ai,Vi)==eye(5) and m(m(sel(Jm,I),Vi),A)==Jm
Z=m(Jm,tr(M));beta=max(F(Z[i][i]+Z[j][j]-2*Z[i][j],96**2)for i,j in combinations(range(6),2))
assert beta==F(7,16)
assert all(any((2*v-w)%96 for v,w in zip(S[1],x))for x in S)
k=read('contact.txt');assert len(k)==1 and len(k[0])==64
assert m(m(k,J),tr(k))==[[12]]
assert m(k,tr(A))==[[0,0,0,-1,0]]
d=sub([M[0]],[M[4]])
assert F(m(m(d,J),tr(d))[0][0],96**2)==beta
delta=F(36)*F(185,64)**59/(32768*3**28)
v=json.loads((P/'data/results.json').read_text())
assert F(v['center_density_squared_floor'])==delta>F(3**33*401**59,2**296*25*13**61)
print('PASS literal period determinant, affine shifts, exact gain ratio, periodic not lattice')
