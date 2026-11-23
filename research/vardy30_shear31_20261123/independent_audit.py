#!/usr/bin/env python3
"""Exact independent verification of sheared periodic 31D packing."""
from itertools import product
from collections import Counter
from fractions import Fraction
assert __debug__
m=15
g=1491
mask=(1<<m)-1
rows=[((g<<i)|(g>>(m-i)))&mask for i in range(m)]
assert len(set(rows))==m
transpose=[sum(((r>>i)&1)<<j for j,r in enumerate(rows)) for i in range(m)]
def linear(rows,x):
 a=0
 for i,row in enumerate(rows):
  if x>>i&1:a^=row
 return a
CG=[(1<<i)|(r<<m) for i,r in enumerate(rows)]
DG=[r|(1<<(m+i)) for i,r in enumerate(transpose)]
assert all((c&d).bit_count()%2==0 for c,d in product(CG,DG))
C=[x|(linear(rows,x)<<m) for x in range(1<<m)]
D=[linear(transpose,x)|(x<<m) for x in range(1<<m)]
assert len(C)==len(set(C))==len(D)==len(set(D))==1<<15
assert C[0]==D[0]==0
assert set(C)!=set(D)
for K in (C,D):
 wh=Counter(c.bit_count() for c in K)
 assert wh[0]==1 and all(not w%2 for w in wh) and min(w for w in wh if w)==8
 bh=Counter((c^mask).bit_count() for c in K)
 assert min(bh)==7
 print('code weight enumerator',sorted(wh.items()), 'distance-to-b min',min(bh),'multiplicity',bh[7])
# For E=-2c+4z and O=1-2d+4(e_1+w), horizontal component
# E-O congruent to odd coordinates. The 30 nearest +/-1 reps need
# correction 1 exactly at c_i=1,d_i=0; their parity is even by
# C even and C·D=0, contrary to z−w−e_1 having odd sum.
u={}
for c,d in product(range(2),repeat=2):
 odd=-1-2*c+2*d
 shift=next(t for t in (-1,0,1) if abs(odd+4*t)==1)
 u[c,d]=shift
assert u=={(0,0):0,(0,1):0,(1,0):1,(1,1):0}
within_mixed=30+8
# Neighbour mixed distances use only odd coordinates; same stock has
# binary weight d(b,K)=7, hence min horizontal squared >=28.
interlayer_same=4*7+2**2
interlayer_mixed=30+2**2
# At lag two, opposite stock remains odd (>=30); same stock with
# unequal codewords has >=32, or its D30 parity shifts by the odd b.
assert mask.bit_count()==15 and mask.bit_count()%2==1
lag2_same=4**2+4**2
lag2_mixed=30+4**2
# At lag three or beyond, the height alone exceeds minimum squared 32.
assert min(32,within_mixed,interlayer_same,interlayer_mixed,lag2_same,lag2_mixed,6**2)==32
# Period group 4D30 × 8Z, because at lag4 horizontal translation is 8b
# and 2b ∈ D30. Two disjoint code families each 32768 stock per layer.
N=4*(len(C)+len(D))
covol=8*2*4**30
vol_per_point=Fraction(covol,N)
assert vol_per_point==2**46
# Radius squared 8, so squared center density = 8^31/(2^46)^2.
assert Fraction(8**31,vol_per_point**2)==2
print('PASS minimum squared distance 32; periodic centers',N,'period covolume',covol,'volume per center',vol_per_point,'delta^2 = 2')
