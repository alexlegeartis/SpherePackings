import runpy,json
from pathlib import Path
from itertools import combinations, product
from math import isqrt
from fractions import Fraction as Q
import numpy as np
from sympy import Matrix
from hashlib import sha256
import sys
P=Path(__file__).resolve().parent
H=P.parent/'gamma72_depth6_20261122'
pins={'verify.py':'ca65ade908d1399acb93c549f304dcf565f731d347a9ebe525cf1e5a28c7d2c6',
      'certificate.json':'60d65708796ece1a76afd974d77e509fbe01dd2ceef37fbffbf6dd7dfffd8870',
      'numerators.json':'c6a30964135f477fc361bd2c9e6350a4041524d763db2fe901fc2e8abde5fd59',
      'results.json':'46393d59baa87284c98bb4f45920f401f78a0ad66c67c6e612739d397654ba00'}
for f,h in pins.items():
 assert sha256((H/f).read_bytes()).hexdigest()==h, f
d=runpy.run_path(str(H/'verify.py'))
phases=json.loads((P/'phases.json').read_text())
W0=np.array(phases['old_D4']+phases['new'],dtype='i8')
assert W0.shape==(7,72) and np.array_equal(W0,W0%2)
Base=d['V0']; U=Base[:,:3]/2
V0=Base.row_join(U).col_join(U.T.row_join(Base[:3,:3]))
def depth(raw,need):
 s,t,b=d['syndrome'](raw); o=d['orth']; M=d['M4']
 if not(s or t or b):return 0
 if need<=2:return 2
 if not(s and t and s!=t):
  assert need<=4,(s,t,b,need)
  if not(s or t):
   assert b
   return 4
  a=s or t
  if b in M[a]:return 2
  if not o(a,b):return 3
  return 4
 if need<=4:return 4
 tri=(s,t,s^t)
 if any(b==x^y for u,w in combinations(tri,2) for x in M[u] for y in M[w]):return 4
 if need<=5:return 5
 if any(not o(k,b^j) for a in tri for k in tri if a!=k for j in M[a]):return 5
 return 6
rows={x['dimension']:x for x in json.loads((P/'results.json').read_text())['rows']}
assert set(rows)=={77,78,79}
for n,expdet in ((77,6),(78,Q(27,4)),(79,Q(27,4))):
 V=V0[:n-72,:n-72]; W=W0[:n-72]
 assert V.is_positive_definite and V.det()==expdet
 v2=np.array((2*V).tolist(),dtype='i8'); assert Matrix(v2.tolist())==2*V
 bounds=[isqrt(int(8*V.inv()[i,i])) for i in range(n-72)]
 seen=set()
 for z in product(*(range(-b,b+1) for b in bounds)):
  if not any(z):continue
  a=np.array(z,dtype='i8'); v=int(a@v2@a)
  if v>=16:continue
  raw=a@W%2; need=(17-v)//2
  dep=depth(raw,need); assert dep>=need, (n,z,v,dep,need)
  seen.add((tuple(a%2),v))
 assert len(seen)=={77:67,78:147,79:295}[n]
 N=2*Matrix(W.tolist())*d['Ti']
 assert all(k.q==1 for k in N)
 G=d['GG']; new=G.row_join(G*N.T/2).col_join((N*G/2).row_join(N*G*N.T/4+V))
 assert new.det(method='domain-ge')==expdet
 print('PASS full coefficient box',n,'bounds',bounds,'distinct',len(seen),'det',expdet,flush=True)
 pi=d['pi']; ball=(pi**n/d['factorial'](n//2)**2 if n%2==0 else 2**(n+1)*pi**(n-1)/d['prod'](range(1,n+1,2))**2)
 gt=ball*2**n/expdet
 row=rows[n]
 assert str(expdet)==row['gram_determinant']
 assert Q(row['center_density_squared'])==Q(2**n,expdet)
 assert Q(row['occupied_density_floor'])**2<gt
 for ref,gain in row['comparisons'].values():
  g=Q(gain)
  assert (1+g/100)**2<=Q(ref)/expdet<(1+(g+Q(1,10**6))/100)**2
 M=4*new
 assert all(z.q==1 for z in M)
 gram={'dimension':n,'minimum_squared':8,'denominator':4,
       'numerator':[[int(z) for z in r] for r in M.tolist()]}
 dest=P/f'gram{n}.json'
 if '--emit' in sys.argv:
  dest.write_text(json.dumps(gram,separators=(',',':'))+'\n')
 else:
  assert gram==json.loads(dest.read_text())
 print('PASS exact Gram, density floor and comparisons',n,flush=True)
