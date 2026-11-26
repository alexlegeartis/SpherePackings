#!/usr/bin/env python3
"""MW92 95D conditional tail construction. Exact rational arithmetic."""
from fractions import Fraction as F
from itertools import product,permutations
from math import isqrt,prod
from decimal import Decimal
import json
assert __debug__
A=F('11.00098526'); p0=16-A; t=16-3*p0
D={6:F('10.58'),7:F('10.545'),9:F('10.487')}
assert all(k*d**46<A**46 for k,d in D.items())
N=lambda x,y:x*x+x*y+y*y
near=[(0,0),(1,0),(0,1),(-1,1),(1,-1),(1,1)]
far=[(0,-1),(-1,0),(-1,2),(0,2),(2,-1),(2,0)]
S=near+far
def diff(c,e):return N(c[0]-e[0],c[1]-e[1])
assert all(N(3*x-1,3*y-1) in (3,12) for x,y in near)
assert all(N(3*x-1,3*y-1)==21 for x,y in far)
assert all(any(diff(c,e)==4 for e in near) for c in far)
assert {diff(c,e) for c in far for e in far}=={0,1,4,7,9}
assert min(D.values())>9 and A>11 and t>1
def residues(fn,n,loss):
 rep={}
 for c in S:rep.setdefault(fn(*c),c)
 assert len(rep)==n and set(near)<=set(rep.values())
 for c in S:assert diff(c,rep[fn(*c)]) in (0,loss)
for w in (3,5):
 assert (w*w-w+1)%7==0
 residues(lambda x,y:(x+w*y)%7,7,7)
residues(lambda x,y:(x%3,y%3),9,9)
def orbit(x,y):
 out=set()
 for _ in range(6):out.add((x,y));x,y=-y,x+y
 return out
O=[orbit(*z) for z in [(1,0),(1,1),(2,0),(2,1),(1,2),(3,0)]]
for n,ids in [(1,[0]),(3,[1]),(4,[2]),(7,[3,4]),(9,[5])]:
 # N>=3x²/4,3y²/4, hence |x|,|y|<=5 suffice here.
 assert {(x,y) for x in range(-5,6) for y in range(-5,6) if N(x,y)==n}==set.union(*(O[i] for i in ids))
for w in (3,5):
 for oc in (3,4):
  assert len({(x+w*y)%7==0 for x,y in O[oc]})==1
 assert sum(all((x+w*y)%7==0 for x,y in O[oc]) for oc in (3,4))==1
def det(g):
 if not g:return F(1)
 return sum((-1)**j*z*det([r[:j]+r[j+1:] for r in g[1:]]) for j,z in enumerate(g[0]))
def short(g):
 n=len(g);d=det(g)
 assert all(det([r[:j] for r in g[:j]])>0 for j in range(1,n+1))
 b=[isqrt(int(16*det([r[:i]+r[i+1:] for r in g[:i]+g[i+1:]])/d)) for i in range(n)]
 for z in product(*(range(-x,x+1) for x in b)):
  if not any(z) or next(x for x in z[::-1] if x)<0:continue
  q=sum(z[i]*g[i][j]*z[j] for i in range(n) for j in range(n))
  if q<16: yield z,q
def generic(k):
 d=D[k];p=F(5);u=d/4 if k==6 else p/2;c=16-p-2*u if k==6 else 16-d
 return [[p,p/2,-u],[p/2,p,-u],[-u,-u,c]]
def check(k):
 g=generic(k);cnt=0
 for (i,j,b),q in short(g):
  if b==0:h=A if N(i,j)==1 else t if N(i,j)==3 else F(0)
  elif b==1:h=D[k] if (i,j) in near else F(1) if (i,j) in far and (k!=6 or (i,j)!=far[0]) else F(0)
  else:h=F(0)
  assert q+h>=16,(k,i,j,b,q,h)
  cnt+=1
 return det(g),cnt
V=[(0,0),(1,0),(0,1)];moved=set()
for v0,v1,v2 in permutations(V):
 def act(c):return tuple(v0[i]+c[0]*(v1[i]-v0[i])+c[1]*(v2[i]-v0[i]) for i in range(2))
 assert set(map(act,near))==set(near) and set(map(act,far))==set(far)
 moved.add(act(far[0]))
assert moved==set(far)
gb={k:check(k) for k in D}
B0=F(14,3);d3=F('10.74136270');d8=F('10.51475582')
diag=[B0,16-d3,16-d8];g=[[diag[min(i,j)]/(1 if i==j else 2) for j in range(3)] for i in range(3)]
seen={1:set(),2:set()}
for z,q in short(g):
 j=max(i for i,x in enumerate(z) if x)
 if not j:assert z[0]%3 and q+F(34,3)>=16
 else:
  assert z[j]==1 and q+{1:d3,2:d8}[j]>=16
  seen[j].add((z[0]%3,*z[1:j]))
assert len(seen[1])<=3 and len(seen[2])<=8
worst=det(g);prev=F(73305736561980600794262547,10**24)
assert worst==F(652349944595359843,9*10**15) and all(d<worst for d,n in gb.values()) and worst<prev
assert F(1005657,10**6)**2<prev/worst<F(1005658,10**6)**2
C2=F(4)**95/(47**14*worst)
def atan(m,n):return sum(F((-1)**i,(2*i+1)*m**(2*i+1)) for i in range(n))
pi_lo=16*atan(5,80)-4*atan(239,81)
V2=F(2)**96*pi_lo**94/prod(range(1,96,2))**2
q=V2*C2*10**72
lower=isqrt(q.numerator//q.denominator)
print(json.dumps(dict(old_det=str(prev),generic_dets={k:str(v[0])for k,v in gb.items()},short_representatives={k:v[1]for k,v in gb.items()},worst_det=str(worst),squared_center_density=str(C2),occupied_density_lower=str(Decimal(lower)*Decimal(10)**-36),gain_percent_interval=['0.5657','0.5658']),indent=2))
