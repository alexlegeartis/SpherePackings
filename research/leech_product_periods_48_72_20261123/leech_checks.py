from fractions import Fraction as F
from math import prod
from itertools import product
N=196560
# Full signed-shell dot counts for a least rep of norm 4,6,8.
rows=[(93150,94208,9200,0,2),
      (75900,97200,22356,1104,0),
      (65780,94208,32384,4096,92)]
for j,cs in enumerate(rows,2):
 assert sum(cs)==N
 for h in range(1,5 if j==4 else 4):
  rhs=F(N*(8*j)**h*prod(range(1,2*h,2)),prod(range(24,24+2*h,2)))
  assert sum(v*i**(2*h) for i,v in enumerate(cs))==rhs,(j,h)
eig=[N//2-sum(cs[1::2]) for cs in rows]
assert eig==[4072,-24,-24]
assert F(2**24*24,98280+24)==4096
assert F(4096,2**24)*F(8,6)**24<1
assert 2**12*F(6,8)**36<1
assert 2**12*F(4,8)**24<1
for s in range(601,801):
 mu=F(s,100); assert mu**2<18*(mu-4) 
print('PASS exact shell moment system, Hoffman 4096, and comparisons')
