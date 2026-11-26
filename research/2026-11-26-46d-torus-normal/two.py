from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from sympy import Matrix,Rational
from math import sqrt
import json,sys
from functools import cache
def go(path):
 o=json.load(open(path));A,B,g=(o[k] for k in ['A','B','g']);D=A*B-g*g
 aff=[list(map(F,r)) for r in o['aff']];q=list(map(F,o['minima']));n=len(q);d=len(aff[0])-1
 def fl(a):return list(map(float,a))
 C=np.asarray([fl(r[1:]) for r in aff]);e=np.asarray([float(r[0]) for r in aff])
 cert=[]
 def lb(row):
  rec={};cert.append(rec)
  if not d or not any(row[1:]):return row[0]
  r=linprog(fl(row[1:]),A_ub=-C,b_ub=e,bounds=[(None,None)]*d,method='highs')
  if not r.success:return F(0)
  ids=[i for i,b in enumerate(r.ineqlin.marginals) if b<-1e-11]
  if not ids:return F(0)
  mat=Matrix([aff[i][1:] for i in ids]).T
  try:y,par=mat.gauss_jordan_solve(Matrix(row[1:]))
  except ValueError:return F(0)
  if len(par):return F(0)
  y=[F(v) for v in y];assert all(v>=0 for v in y)
  assert all(sum(yj*aff[i][k] for yj,i in zip(y,ids))==row[k] for k in range(1,d+1))
  rec.update({str(i):str(yj) for yj,i in zip(y,ids) if yj});return row[0]-sum(yj*aff[i][0] for yj,i in zip(y,ids))
 best=F(0);cand=[]; res=[]
 for i,h in enumerate(q):
  if not h:
   z=lb(aff[i]);assert z>0; res.append((str(z),'0','0'));continue
  a=lb(aff[i]);b=F(0)
  if a>0:cost=h
  else:
   b=lb([x+y for x,y in zip(aff[i],aff[n+i])])
   cost=max(F(0),h-2) if b>0 else F(0)
  val=F(4,D)*((6-cost)/4)**46;best=max(best,val)
  if val>F(1,3)*F(13,9)**46:cand.append((i,o['coords'][i],str(h),str(a),str(b)))
  res.append((str(a),str(b),str(cost)))
 print('2marks',A,B,g,'best <=',sqrt(best),'unexcluded',cand,flush=True)
 with open('lp_'+path,'w') as f:json.dump({'best2':str(best),'rows':res,'survivors':cand,'duals':cert},f)
if __name__=='__main__':
 for path in sys.argv[1:]:go(path)
