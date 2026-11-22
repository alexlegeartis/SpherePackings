import runpy,json
from pathlib import Path
from itertools import combinations, product
from math import isqrt
from fractions import Fraction as Q
import numpy as np
from sympy import Matrix
HERE=Path(__file__).resolve().parent
H=HERE.parent/'gamma72_depth6_20261122'; d=runpy.run_path(str(H/'verify.py'))
phases=json.loads((HERE/'phases.json').read_text())
W=np.vstack([d['rawD4']]+[p['raw'] for p in phases])
V=Matrix([[4,-2,0,0,1],[-2,4,-2,-2,0],[0,-2,4,0,0],[0,-2,0,4,0],[1,0,0,0,4]])/2
V=V.row_join(Matrix([1,0,0,0,-1])/2).col_join(Matrix([[1,0,0,0,-1,4]])/2)
v2=np.array((2*V).tolist(),dtype='i8')
Vinv=V.inv()
assert V.is_positive_definite and V.det()==10
bounds=[isqrt(int(8*Vinv[i,i])) for i in range(6)]
seen=set()
for z in product(*(range(-b,b+1) for b in bounds)):
 if not any(z):continue
 a=np.array(z,dtype='i8')
 v=int(a@v2@a)
 if v>=16:continue
 raw=a@W%2
 s,t,b=d['syndrome'](raw)
 if not any(z[4:]):
  if v==4:d['prove'](raw,6)
  if v==8:
   if s==0 and t==0:d['prove'](raw,4)
   elif s==t or not(s and t):assert b not in d['M4'][s or t] and d['orth'](s or t,b)
  if v==12:assert s or t or b
  assert v in (4,8,12)
 else:
  assert s and t and s!=t
  assert all(b!=x^y for u,w in combinations((s,t,s^t),2) for x in d['M4'][u] for y in d['M4'][w])
  assert v>=6 or tuple(z[:4])==(0,0,0,0) and v==4 and sum(map(abs,z[4:]))==1
  if v==4:d['prove'](raw,6)
  seen.add((tuple(a%2),v))
assert len(seen)==105
N=2*Matrix(W.tolist())*d['Ti']
assert all(k.q==1 for k in N)
G=d['GG']; new=G.row_join(G*N.T/2).col_join((N*G/2).row_join(N*G*N.T/4+V))
assert new.det(method='domain-ge')==10
pi=d['pi']; ball=pi**78/d['factorial'](39)**2
gt=ball*2**78/10
from decimal import Decimal, localcontext
with localcontext() as ctx:
 ctx.prec=28
 x=(Decimal(gt.numerator)/Decimal(gt.denominator)).sqrt()
 print('78D exact preflight:',len(seen),'tilted parity/norm cases, exact min 8, det 10, density',x,'gain',(Decimal(12)/10).sqrt()*100-100)
for row in json.loads((HERE/'results.json').read_text())['rows']:
 n=row['dimension']
 k=row['gram_determinant']
 assert n in (77,78)
 gg=new[:n,:n]
 if n==77:
  # The earlier generator is the upper-left principal submatrix.
  assert phases[0]['raw']==json.loads((HERE/'phase.json').read_text())['raw']
 assert gg.det(method='domain-ge')==k
 assert Matrix(json.loads((HERE/f'gram{n}.json').read_text())['gram_numerator'])==4*gg
 v=Q(2**n,k)
 assert v==Q(row['center_density_squared'])
 ball=(pi**n/d['factorial'](n//2)**2 if n%2==0 else
       2**(n+1)*pi**(n-1)/d['prod'](range(1,n+1,2))**2)
 assert Q(row['occupied_density_floor'])**2<ball*v
 for oldkey,gainkey in (
  ('preceding_repository_gram_determinant','gain_over_preceding_repository_percent_floor'),
  ('named_external_gram_determinant','gain_over_named_external_percent_floor')):
  g=Q(row[gainkey])
  assert (1+g/100)**2<=Q(row[oldkey],k)<(1+(g+Q(1,10**6))/100)**2
 print(f'PASS exact density floor and comparisons, {n}D')
