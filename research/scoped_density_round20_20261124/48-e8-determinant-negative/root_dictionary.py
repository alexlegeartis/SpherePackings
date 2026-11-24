#!/usr/bin/env python3
import json
import subprocess
from pathlib import Path
from struct import unpack_from
from fractions import Fraction as F
from tempfile import TemporaryDirectory
from collections import Counter
ROOT=Path(__file__).resolve().parent
TARGET=F(38246144,2401)
HERM=F(389,2000)
B6=38246144
N0=2**48
def check_lifts():
 from math import prod
 def q(x):return (x&(x>>4)).bit_count()%2
 def b(x,y):return q(x^y)^q(x)^q(y)
 def dot(x,y):return sum(a*b for a,b in zip(x,y))
 R=[tuple(map(int,l.split())) for l in (ROOT/'roots.in').read_text().splitlines()]
 A=[x for x in range(256) if q(x)]
 assert len(R)==256 and len(A)==120 and all(len(r)==8 for r in R)
 for x,v in enumerate(R):
  if x not in A:
   assert not any(v)
   continue
  assert dot(v,v)==8 and len({k%2 for k in v})==1 and sum(v)%4==0
  for y in A:
   z=dot(v,R[y])
   assert z%4==0 and (z//4)%2==b(x,y)
 span={0}
 for x in A:
  span|={x^u for u in span}
 assert len(span)==256
 base=[1,16,2,32,4,64,8,128]
 order=[135,64,35,16,9,4,2,1]
 extra=[21,25,29,17]
 assert all(x in A for x in extra)
 def tr(x,r):return x^r if b(x,r) else x
 def comp(x):
  for r in extra:x=tr(x,r)
  return x
 assert all(comp(j)==j for j in base[:5])
 for i, x in enumerate(base):
  ts=[r for r in A if all(not b(r,t) for t in base[:i])]
  orbit={x}; stack=[x]
  for u in stack:
   w=[tr(u,r) for r in ts]
   if i==5:w.append(comp(u))
   for v in w:
    if v not in orbit:orbit.add(v);stack.append(v)
  assert len(orbit)==order[i]
 size=prod(order)
 assert size==prod(((1<<(2*m-1))+(1<<(m-1))-1)*(1<<(2*m-2)) for m in range(1,5))==348364800
 print('E8 root dictionary and surjectivity of liftable orthogonal group PASS',flush=True)

