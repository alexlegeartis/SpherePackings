#!/usr/bin/env python3
"""Independent reverse-LDL audit of the finite order-three branch."""
from fractions import Fraction as F
from math import prod
import json
from pathlib import Path

def ldl(A):
 n=len(A);D=[];L=[[F(i==j)for j in range(n)]for i in range(n)]
 for i in range(n):
  for j in range(i):
   L[i][j]=(A[i][j]-sum(L[i][k]*L[j][k]*D[k]for k in range(j)))/D[j]
  D.append(A[i][i]-sum(L[i][k]**2*D[k]for k in range(i)))
  assert D[i]>0
 return L,D
def enumerate_short(G):
 L,D=ldl(G);n=len(D)
 v=[0]*n
 def layer(j,room):
  if j<0:
   if any(v):yield tuple(v),16-room
   return
  t=sum(L[i][j]*v[i]for i in range(j+1,n))
  m=(-t).numerator//(-t).denominator
  k=m
  while D[j]*(k+t)**2 <room:
   v[j]=k;yield from layer(j-1,room-D[j]*(k+t)**2);k-=1
  k=m+1
  while D[j]*(k+t)**2 <room:
   v[j]=k;yield from layer(j-1,room-D[j]*(k+t)**2);k+=1
 yield from layer(n-1,F(16))

here=Path(__file__).resolve().parent
out=json.loads((here/'results.json').read_text())
G=[[F(q)for q in row]for row in out['branch_B_gram']]
L,D=ldl(G)
d=list(map(F,out['branch_B_depths']))
c=[set()for _ in d]
r=F('11.000980284')
n=0
for z,q in enumerate_short(G):
    j=max(i for i,a in enumerate(z)if a)
    if z[j]<0:continue
    n+=1
    assert z[j]==1 and q+d[j]>=16
    c[j].add((z[0]%3,)+z[1:j] if j else ())
assert n==out['branch_B_short_pairs']==31
assert list(map(len,c))==out['branch_B_counts']==[1,3,7,17]
for j in range(1,4):
    assert len(c[j])*d[j]**46<r**46
assert prod(D)==F(out['branch_B_det'])
print('PASS exact reverse-LDL:',n,list(map(len,c)),prod(D))
