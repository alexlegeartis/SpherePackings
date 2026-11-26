from fractions import Fraction as F
from itertools import product
from math import sqrt
import json
from pathlib import Path
def bron(N):
 # All graphs are vertex transitive, so fix vertex zero. Maximal-clique recursion,
 # independent of upper-bounding sequential greedy-color search in analyse.py.
 best=1
 def visit(sz,P,X):
  nonlocal best
  if not P: best=max(sz,best);return
  if sz+len(P)<=best:return
  pivot=max(P|X,key=lambda j:len(P&N[j]))
  for v in list(P-N[pivot]):
   visit(sz+1,P&N[v],X&N[v]);P.remove(v);X.add(v)
 visit(1,set(N[0]),set());return best
def check(path):
 o=json.load(open(path)); hh=json.load(open('lp_'+path))
 D=o['A']*o['B']-o['g']**2;K=list(map(tuple,o['keys']));n=len(K)
 aff=[list(map(F,r)) for r in o['aff']];d=len(aff[0])-1
 data=iter(hh['duals']); costs={}
 def verify(row,claim):
  dual=next(data); mult={int(k):F(t) for k,t in dual.items()}
  assert all(v>=0 for v in mult.values())
  if mult or not any(row[1:]):
   assert all(sum(v*aff[k][j] for k,v in mult.items())==row[j] for j in range(1,d+1))
   assert row[0]-sum(v*aff[k][0] for k,v in mult.items())==claim
  else: assert claim<=0 # no dual needed; each population is nonnegative
 for i,(key,(a,b,c),h) in enumerate(zip(K,hh['rows'],map(F,o['minima']))):
  a,b,c=map(F,(a,b,c));verify(aff[i],a)
  if not h:assert a>0
  if h and a<=0:verify([x+y for x,y in zip(aff[i],aff[n+i])],b)
  assert c==(h if a>0 else max(F(0),h-2) if b>0 else F(0))
  costs[key]=costs[tuple((-x)%D for x in key)]=c
 assert next(data,None) is None
 G=sorted(costs);assert len(G)==D
 best=(F(0),)
 for b in sorted(set(costs.values())):
  if b>=6:continue
  sc=((6-b)/4)**46 / D
  if D**2*sc<=best[0]:continue
  adj=[{j for j,(u,v) in enumerate(G) if i!=j and costs[(u-x)%D,(v-y)%D]<=b} for i,(x,y) in enumerate(G)]
  s=bron(adj);best=max(best,(s*s*sc,s,b))
 assert best[0]<=F(13,9)**46/3
 return {'gram':[o['A'],o['g'],o['B']],'det':D,'size':best[1],'beta':str(best[2]),'delta_squared_upper':str(best[0]),'delta_approx':sqrt(best[0])}
if __name__=='__main__':
 cases=[(a,b,g) for a in range(6,11,2) for b in range(a,13,2) for g in range(a//2+1) if a*b-g*g<=68]
 assert len(cases)==20
 assert 68<F(12)*(F(27,26)**46)<69
 res=[check(f'mixed_{a}_{b}_{g}.json') for a,b,g in cases]
 assert sum(F(r['delta_squared_upper'])==F(13,9)**46/3 for r in res)==1
 print(len(res),'reduced Gram types audited; all <=CS46; 6,3,6 uniquely equals baseline')
 with open('summary46.json','w') as f:json.dump(res,f,indent=1)
