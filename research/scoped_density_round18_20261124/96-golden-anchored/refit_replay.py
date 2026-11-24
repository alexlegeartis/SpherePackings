"""Exact all-degree refit against the A=0 genus-two dot budgets."""
from pathlib import Path
import sys
# Reuse the unchanged round-17 interval arithmetic and span checker.
prior=Path(__file__).resolve().parents[2]/'scoped_density_round17_20261124'/'96-product-two-point'
sys.path.insert(0,str(prior))
from replay import *

N=1218201435360000
src=json.loads(Path(__file__).with_name('refit-exact.json').read_text())
D={(p,b):F(w) for p,b,w in src}
assert len(src)==len(D)
assert D.get((8,0),0)==0
assert 2+sum(D.values())==N
assert 288+sum(w*p*p for (p,b),w in D.items())==F(3,2)*N
assert sum(w*b*b for (p,b),w in D.items())==F(15,2)*N
for (p,b),w in D.items():
 for q,s in span_ok(p,b):
  assert w/((2 if p else 1)*(2 if b else 1))<=size(0,q,s)
table=[]; tail=iv(0)
for (p,b),w in D.items():
 assert w>0
 u,v=angles(p,b)
 assert all(max(abs(z) for z in x)<S for x in (u,v))
 table.append((iv(w/2),normed(u,K),normed(v,K)))
 tail=add(tail,mul(iv(w/2),add(envelope(u),envelope(v))))
assert tail[1]<S//400
lb=2*S
for i in range(K):
 for j in range(i+1):
  if (i+j)%2 or (i,j) in ((0,0),(1,1),(2,0)):continue
  x=iv(2)
  for w,r,t in table:x=add(x,mul(w,add(mul(r[i],t[j]),mul(r[j],t[i]))))
  assert x[0]>S//5,(i,j,float(F(x[0],S)))
  lb=min(lb,x[0])
print('exact finite minimum',float(F(lb,S)),'tail',float(F(tail[1],S)))
import jacobi_replay as j
ta=j.theta(0); m,_=j.mj(5,5,ta)
c,_=j.mj(6,5,ta,{5:{r:m(6,r)[0] for r in range(-12,13)}})
for p in range(8):
 used=F(N)*sum(w for (r,b),w in D.items() if r==p)/(2 if p else 1)
 assert used<c(6,p)[0]
 print('absolute dot',p,', used fraction',float(used/c(6,p)[0]))
