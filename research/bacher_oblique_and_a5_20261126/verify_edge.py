#!/usr/bin/env python3
"""Exact A5-to-D5 edge exclusion and A5 map count. Standard library. Read verify_glue.py for geometry."""
from itertools import product
from fractions import Fraction as F
W=[0,3,4,3,2,3,2,7,8,7,4,7,2,7,2,3,4,3,8,3,2,7,2,3,4,7,8,7,2,3,2,7,4,3,8,3,2,7,2,3,4,7,8,7,2,3,2,7,4,7,8,7,6,7,6,3,4,3,8,3,6,3,6,7]
H=[0,2,8,10,16,18,24,26,32,34,40,42,48,50,56,58]
def add(a,b):return ((a^b)^(((a&b)&341)<<1))&1023
T=[a for a in range(1024) if not add(a,a)]
roots={x:[z for z in range(1024) if add(z,z)==x] for x in T}
def count(C):
 def check(h,off):return all(C[x]>=W[i+off] for x,i in zip(h,H))
 c=n=0
 for a,b in product(T,repeat=2):
  v=[0,a,b,add(a,b)]
  if not check(v,0):continue
  for e in T:
   s=v+[add(x,e) for x in v]
   if not check(s,0):continue
   for f in T:
    h=s+[add(x,f) for x in s]
    if not check(h,0):continue
    c+=1
    A=[p for p in roots[a] if check([add(p,x) for x in h],1)]
    B=[p for p in roots[b] if check([add(p,x) for x in h],4)]
    for p,t in product(A,B):n+=check([add(add(p,t),x) for x in h],5)
 return c,n
# On sum(v)=0, v0*v1 <= (sum v_i^2)/3: write common sum s=v0+v1, use the four other coords.
# Hence Q_t=2q-t v0v1 is >=(2/3)q throughout 2<=t<=4. Only q<=48 needed.
V={}
for v in product(range(-6,7),repeat=5):
 q=sum(x*x for x in v)+sum(v)**2
 if q<=48:
  k=(sum((x%4)<<(2*i) for i,x in enumerate(v)),v[0]*v[1])
  V[k]=min(V.get(k,99),q)
C=[99]*1024
for (i,p),q in V.items():C[i]=min(C[i],q//2)
assert count(C)==(720,2880)
# max permitted uniform scaling is (3/2)^(1/5) < 217/200.
assert 2*217**5>3*200**5
ts={F(2),F(4)}
for (i,p),q in V.items():
 if not p:continue
 for w in set(W)-{0}:
  t=F(217*2*q-800*w,217*p)
  if 2<=t<=4:ts.add(t)
T0=sorted(ts);assert len(T0)==42
# Thresholds only depend on crossing events; pick one rational strictly between successive events.
test=T0+[(a+b)/2 for a,b in zip(T0,T0[1:])]
seen=set()
for t in test:
 a,b=t.as_integer_ratio();C=[99]*1024
 for (i,p),q in V.items():C[i]=min(C[i],(217*(2*q*b-p*a))//(800*b))
 # Identify profiles by which of the actually required thresholds are met.
 sig=bytes(sum(c>=u for u in [2,3,4,6,7,8]) for c in C)
 if sig not in seen:assert count(C)[1]==0;seen.add(sig)
assert len(seen)==19
print("A5: 720 half assignments, 2880 valid maps. Winning edge: 42 breakpoints, 83 representatives, 19 profiles, zero maps.")
