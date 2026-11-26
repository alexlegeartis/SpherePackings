#!/usr/bin/env python3
"""Exact inherited-radius section check. Requires sympy for matrix arithmetic; run from any subdir in the repository."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,combinations
from math import isqrt,prod
from sympy import Matrix
P=next((p/"research/bacher_slip_29_30_20261122" for p in [Path.cwd(),*Path.cwd().parents] if (p/"research/bacher_slip_29_30_20261122/gram.txt").exists()))
import json
G=Matrix([list(map(int,r.split())) for r in (P/"gram.txt").read_text().splitlines() if r.strip()])
w=json.loads((P/"witness.json").read_text());I=[i for i in range(28) if i not in (4,27)];N=[4,27];Gi=G.inv()
assert list(-6*Gi[:,4])==w["direction"] and Gi.extract(N,N)==Matrix([[2,-1],[-1,2]])/3
assert all(not(m>>i&1)for m in w['half_labels'] for i in N)
S=G.extract(I,I);assert S.det()==3*2**26
U=Matrix([[0,0,0,0,0,0,0,1,0,0,-1,0,0,0,-1,1,1,-1,-1,0,-1,0,0,0,1,0],[0,0,0,0,-1,0,0,-1,0,-1,1,0,0,1,0,-1,0,0,0,0,0,0,-1,0,-1,0],[0,0,0,0,0,0,1,0,0,1,1,-1,-1,-1,0,0,0,0,1,1,0,0,0,0,0,0],[0,0,0,-1,-1,0,1,0,0,0,-1,-1,0,-1,-1,0,1,0,0,1,-1,-1,0,-1,0,0],[0,0,0,-1,-1,0,1,0,0,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,-1,0,1],[0,0,0,-1,-1,0,1,0,0,0,-1,-1,0,-1,-1,0,0,0,0,0,0,-1,0,-1,0,1],[0,0,0,-1,0,1,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,1,-1,0,0,1,0],[0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,-1,-1,0],[0,0,0,1,0,-1,0,0,0,1,0,0,-1,0,0,1,0,1,0,0,0,0,0,1,1,-1],[0,0,0,0,0,0,0,-1,0,-1,1,0,1,0,1,-1,0,0,1,0,0,0,0,0,-1,0],[0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,1,0,-1,0,0,0,1,1,0,0,1],[0,0,0,0,1,-1,1,0,0,1,0,0,0,-1,1,1,-1,0,1,0,1,-1,1,0,1,1],[0,0,0,-1,0,0,1,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,-1,-1,1],[0,0,0,0,1,-1,-1,-1,0,1,0,-1,1,-1,0,0,1,1,1,1,0,-1,1,-1,0,0],[0,0,0,-1,0,1,1,0,0,0,0,0,0,-1,0,0,-1,0,0,0,0,1,0,0,-1,0],[0,0,0,0,1,-1,0,-1,0,-1,0,0,0,0,0,0,1,0,0,0,-1,0,0,1,0,0],[0,0,0,1,1,-2,0,-1,0,1,0,-1,0,-1,-1,1,1,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,-1,-1,-1,0,0,0,0,0,-1,1,-1,-1,1,0,0,0,0,0,1,0,-1],[0,0,0,0,0,0,0,0,0,-1,1,0,2,1,1,0,1,0,1,0,1,-1,0,-1,-1,1],[0,0,1,0,1,-1,-1,-1,0,0,0,0,1,0,1,0,0,0,1,0,0,0,1,0,0,0],[0,0,0,0,0,0,-1,0,1,-1,0,1,2,1,0,-1,0,0,-1,-1,0,0,-1,0,0,0],[0,0,1,1,0,0,0,1,1,0,0,1,1,1,0,1,0,-1,0,-1,1,0,0,1,1,1],[0,0,0,-1,0,1,1,1,-1,0,0,0,-1,0,0,1,-1,0,0,-1,0,0,0,0,0,0],[0,0,0,0,-1,0,1,0,0,0,0,0,0,0,-1,0,-1,0,-1,-1,1,1,-1,0,-1,1],[1,1,0,1,1,0,-1,1,0,1,1,0,0,0,1,1,1,1,1,1,0,0,2,0,0,0],[0,1,0,0,1,0,0,0,-1,1,0,0,-1,0,0,1,0,0,0,1,0,0,1,1,1,-1]])
assert abs(U.det())==1
Z=U*(6*S.inv())*U.T;assert all(x.is_Integer for x in Z)
# Integer fraction-free positive Gram decomposition and rigorous exhaustive enumeration:
def shell(Z,cap):
 n=Z.rows; G=[list(map(int,Z.row(i))) for i in range(n)];B=[r[:] for r in G];ds=[];co=[];de=[];p=1
 for i in range(n):
  z=B[i][i];assert z>0;ds.append(z);co.append([B[j][i] for j in range(n)]);de.append(z*p)
  for j in range(i+1,n):
   for k in range(i+1,n):
    a=B[j][k]*z-B[j][i]*B[i][k];assert not a%p;B[j][k]=a//p
  p=z
 scale=1<<20;x=[0]*n;hist={}
 def dfs(i,m):
  if i<0:
   q=sum(x[j]*sum(g*y for g,y in zip(row,x)) for j,row in enumerate(G))
   if q<=cap:hist[q]=hist.get(q,0)+1
   return
  s=sum(co[i][j]*x[j] for j in range(i+1,n));t=isqrt(((m+1)*de[i]-1)//scale)
  for a in range(-((t+s)//ds[i]),(t-s)//ds[i]+1):
   x[i]=a;dfs(i-1,m-(ds[i]*a+s)**2*scale//de[i])
 dfs(n-1,cap*scale);return hist
assert shell(Z,8)=={0:1,8:1638}
# Normal Gram matrices: restrict a dual row to core zero, then lift its two B28 normal coordinates and h vertical coordinates.
rks=[[F(1),F(2,3),F(1,3)],[F(1),F(6,11),F(4,11),F(2,11)],[F(1),F(22,35),F(12,35),F(8,35),F(4,35)]]
mats=[(3,[[2,-1],[-1,2]]),(66,[[48,-24,12],[-24,45,-6],[12,-6,36]]),(210,[[180,-90,60,60],[-90,150,-30,-30],[60,-30,132,48],[60,-30,48,132]])]
def rows_below(M,cap):
 Mi=M.inv(); B=[isqrt(int(cap*Mi[i,i])) for i in range(M.rows)]
 return [(v,q) for a in product(*(range(-b,b+1) for b in B)) if any(a) for v in [Matrix(a)] for q in [(v.T*M*v)[0]] if q<=cap]
for h,((den,A),R) in enumerate(zip(mats,rks)):
 M=Matrix(A)/den; assert M.det()==R[-1] and min(q for _,q in rows_below(M,1))==R[1]
 if M.rows>2:
  assert M.det()*min(q for _,q in rows_below(M.inv(),2))==R[-2]
 if M.rows==4:
  # For a reduced basis of any plane, 3 a*b/4<=det; a>=lambda_1, so if det<=R2 both rows have norm<=4R2/(3lambda1).
  v=rows_below(M,F(4,3)*R[2]/R[1]);d2=min(d for ((a,p),(b,q)) in combinations(v,2) for d in [p*q-(a.T*M*b)[0]**2] if d>0);assert d2==R[2]
 D=[F(4,9),F(16,33),F(64,105)][h]
 gam=[1,1,F(4,3),2,4,8]  # gamma_q^q
 for d,target in [(25,F(1,8)),(26,F(1,3)),(27,F(1,3)),(28,F(4,9))]:
  k=28+h-d
  if k<0:continue
  c=[D/(2**k*R[k-q]*(F(4,3)**q/gam[q])) for q in range(k+1) if k-q<len(R)]
  assert max(c)==target
print("exact core dual shell and inherited-radius section thresholds verified")
