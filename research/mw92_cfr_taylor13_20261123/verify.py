#!/usr/bin/env python3
"""Strict analytic refinements via exact odd Taylor–Dirichlet moments."""
from fractions import Fraction as Q
from itertools import product,permutations
from math import prod,factorial,isqrt
from hashlib import sha256
from pathlib import Path
from decimal import Decimal,localcontext
import json
from moments import cover
assert __debug__
P=Path(__file__).resolve().parent
parent=47**7
d={1:Q('11.00098526'),3:Q('10.74136270'),7:Q('10.54532299'),8:Q('10.51475582'),12:Q('10.42248121')}
def atan(inv,m):return sum(Q((-1)**i,(2*i+1)*inv**(2*i+1)) for i in range(m))
pi_lo=16*atan(5,80)-4*atan(239,81)
pi_hi=16*atan(5,81)-4*atan(239,80)
assert pi_lo<pi_hi
c3,c13=cover(92,3),cover(92,13)
assert c13>c3
for k,a in d.items():assert k*pi_hi**46*a**46<factorial(46)*parent*c13
def det(A):
 n=len(A)
 return sum((-1)**sum(u[i]>u[j]for i in range(n)for j in range(i+1,n))*prod(A[i][j]for i,j in enumerate(u))for u in permutations(range(n)))
def gram(diag):
 return [[diag[min(i,j)]/(1 if i==j else 2) for j in range(len(diag))]for i in range(len(diag))]
def short(G):
 n=len(G)
 assert all(det([row[:i]for row in G[:i]])>0 for i in range(1,n+1))
 D=det(G)
 cof=[det([[x for j,x in enumerate(row)if j!=i]for k,row in enumerate(G) if k!=i]) for i in range(n)]
 B=[isqrt(int(16*a/D))for a in cof]
 for z in product(*(range(-a,a+1)for a in B)):
  if next((v for v in z[::-1] if v),0)<=0:continue
  q=sum(G[i][j]*z[i]*z[j]for i in range(n)for j in range(n))
  if q<16:yield z,q
p=16-d[1]; t=16-3*p
assert d[1]-t/3==Q(32,3)
# Either d(Ta)^2>=t, or a0 has strict depth>32/3 and even/3 norm lattice
# so depth>=34/3. In case B the shift has exact order3.
def pinned_json(name, digest):
 b=(P/name).read_bytes()
 assert sha256(b).hexdigest()==digest
 return json.loads(b)
old0=pinned_json('old.json','3da0691c99ab510a124930a7771972497b0047dc7e62e1f91ab6219e0c78a131')
old1=pinned_json('prior95.json','5c95f72831e8c0292e5693cc2246f71214ed1b969dcc007c8d3d314a304df296')
def branchA(r):
 G=gram([p,p,16-d[12]][:r])
 c=set()
 for z,q in short(G):
  x,y=z[:2]
  i=z[2] if r>=3 else 0
  if not i:
   N=x*x+y*y+x*y
   assert (N==1 and q+d[1]>=16) or (N==3 and q+t>=16)
  else:
   assert i==1 and q+d[12]>=16
   c.add((x,y))
 assert len(c)<=(12 if r==3 else 0)
 return G
def branchB(r):
 G=gram([Q(14,3),16-d[3],16-d[8]][:r])
 seen={1:set(),2:set()}
 for z,q in short(G):
  pos=max(i for i,x in enumerate(z)if x)
  if not pos:
   assert z[0]%3 and q+Q(34,3)>=16
  else:
   assert z[pos]==1 and q+d[{1:3,2:8}[pos]]>=16
   seen[pos].add((z[0]%3,*z[1:pos]))
 assert len(seen[1])<=3 and len(seen[2])<=8
 return G
def ball2(n):
 return pi_lo**n/factorial(n//2)**2 if n%2==0 else 2**(n+1)*pi_lo**(n-1)/prod(range(1,n+1,2))**2
def floor_sqrt(q,digits=16):
 with localcontext() as c:
  c.prec=45
  x=(Decimal(q.numerator)/Decimal(q.denominator)).sqrt()
  a=Decimal(1).scaleb(x.adjusted()+1-digits)
 t=q/Q(a)**2
 k=isqrt(t.numerator//t.denominator)
 assert Q(k*a)**2<=q<Q((k+1)*a)**2
 return f'{a*k:.{digits-1}E}'
def row(n,D,old):
 C=Q(4)**n/(parent**2*D)
 assert C>old
 return {'dimension':n,'det_vertical':str(D),'center_density_squared':str(C),
         'preceding_center_density_squared':str(old),
         'occupied_density_floor':floor_sqrt(C*ball2(n)),
         'gain_percent_approx':100*((float(C/old))**.5-1)}
assert all(q+d[1]>=16 and z==(1,) for z,q in short([[p]]))
ans=[row(93,p,Q(old0['rows'][0]['center_density_squared']))]
for r in (2,3):
 A,B=map(det,(branchA(r),branchB(r)))
 ans.append(row(92+r,max(A,B),Q(old1['rows'][r-2]['center_density_squared'])))
 ans[-1]['both_dets']=[str(A),str(B)]
# 96 general: same nine-offset torus palette, refined d7. D4 in
# physical parity-even coordinates, with an A2 plane in last 3 entries.
z=Q(8,3);s=(3*(16-d[7])-2*z)/4
assert s>=z and 7*z>16
def color(w): return ((w[0]+w[2]+2*w[3])%3,(w[1]+w[2]+w[3])%3)
(u,v),(h,k)=map(color,((1,1,0,0),(1,-1,0,0)))
assert (u*k-v*h)%3
for w in product(range(-2,3),repeat=4):
 if not any(w) or sum(w)%2 or sum(x*x for x in w)>=7:continue
 a,b=color(w)
 q=s*w[0]**2+z*sum(x*x for x in w[1:])+(s-z)*sum(w[1:])**2/3
 threshold=16-d[7] if b else p if a else Q(16)
 assert q>=threshold
ans.append(row(96,4*s*s*z*z,Q(old0['rows'][-1]['center_density_squared'])))
print(json.dumps({'status':'PASS','rows':ans},indent=2))
