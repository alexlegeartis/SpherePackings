"""Replay the added deep-corner and one-third localization bounds using exact leaf weights only."""
from pathlib import Path
import sys,json
from collections import Counter
arg=sys.argv[1:];sys.argv=sys.argv[:1]
import verify as v
from math import sqrt
base=Path(__file__).parent
def check(kind):
 thirds=kind=='thirds'
 states=v.domain() if thirds else {(a,-1) for a in range(1024)}
 reps=[];rem=states.copy()
 while rem:
  a,i=min(rem);reps.append((a,i));g=v.orbit(a,i);assert g<=rem;rem-=g
 file=base/('certificate_'+kind+'.json')
 if arg:
  P=json.loads(Path(arg[0 if thirds else 1]).read_text());assert P['code']==v.C
  file.write_text(json.dumps({'scale':P['scale'],'trees':[[a,i,P['trees'][str(a)][max(0,i)]] for a,i in reps]},separators=(',',':'))+'\n')
 cert=json.loads(file.read_text());S=cert['scale'];assert S==(12288 if thirds else 4096)
 assert [(a,i) for a,i,_ in cert['trees']]==reps
 cnt=Counter()
 def replay(a,l,u,T):
  if T=='interior':
   assert not thirds and a in v.B and sum(z*z for z in u)<S*S
   cnt['inside']+=1;return
  if isinstance(T,dict):
   i=T['j'];assert 0<=i<10 and l[i]<u[i]
   if thirds and 3*l[i]<S<3*u[i]:m=S//3
   else:
    assert (l[i]+u[i])%2==0;m=(l[i]+u[i])//2
   U=u.copy();U[i]=m;L=l.copy();L[i]=m
   replay(a,l,U,T['left']);replay(a,L,u,T['right']);cnt['split']+=1;return
  W={k:int(x) for k,x in T};assert W and min(W)>=0 and max(W)<40 and min(W.values())>0 and len(W)==len(T)
  A=sum(W.values());bound=0
  for i in range(10):
   b=sum(n for k,n in W.items() if ((a^v.C[k])>>i)&1)
   bound+=max((A-b)*y*y+b*(S-y)**2 for y in (l[i],u[i]))
  p,q=(23,9) if thirds else (5,2)
  assert q*bound<=p*A*S*S;cnt['bounded']+=1
 for a,i,T in cert['trees']:
  l=[0]*10;u=[S//2]*10
  if i>=0:l[i]=S//3
  replay(a,l,u,T)
 print(kind,len(states),'states;',len(reps),'roots;',dict(cnt))
for k in ['thirds','radial']:check(k)
for k,d,ref in [(1,.5,9/256),(2,1/(2*sqrt(3)),1/27),(3,1/(4*sqrt(2)),9/256)]:
 b=40/1024*(6/sqrt(13))**k*d
 print(10+k,b,b/ref)

from fractions import Fraction as F
ratios=[(F(40,1024)/ref)**2*F(36,13)**k*norm for k,norm,ref in [(1,F(1,4),F(9,256)),(2,F(1,12),F(1,27)),(3,F(1,32),F(9,256))]]
assert ratios==[F(100,117),F(492075,692224),F(1800,2197)] and all(q<1 for q in ratios)
