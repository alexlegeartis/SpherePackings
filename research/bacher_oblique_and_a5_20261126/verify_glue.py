#!/usr/bin/env python3
"""Check the 27D lattice glue, depending on the published 22D horizontal cost table."""
from itertools import product
from math import prod
from fractions import Fraction as F
W4=[0,3,4,3,2,3,2,7,8,7,4,7,2,7,2,3,4,3,8,3,2,7,2,3,4,7,8,7,2,3,2,7,4,3,8,3,2,7,2,3,4,7,8,7,2,3,2,7,4,7,8,7,6,7,6,3,4,3,8,3,6,3,6,7]
gens=[320,31,2,8] # A5 basis e_i-e_6; base-four little-endian numerator, divide by four.
def add(a,b):return ((a^b)^(((a&b)&341)<<1))&1023
def image(i,gens):
 h=0
 for n,g in zip([i%4,(i>>2)%4,(i>>4)%2,i>>5],gens):
  for _ in range(n):h=add(h,g)
 return h
C=[99]*1024
for v in product(range(-4,5),repeat=5):
 q=sum(a*a for a in v)+sum(v)**2
 if q>32:continue
 j=sum((x%4)<<(2*i) for i,x in enumerate(v));C[j]=min(C[j],q)
# Tail distance squared = C[j]/8. Each smaller candidate has all |v_i|<=4, sufficient through the largest threshold 16.
images=[image(i,gens) for i in range(64)]
assert len(set(images))==64 and all(C[j]>=2*w for j,w in zip(images,W4))
assert sorted(C[j]-2*w for j,w in zip(images,W4))==[0]*57+[8]*6+[16]
# Index 64 over S x sqrt(2)A5 with det S=64 and det(2(I+J))=2^5*6=192:
assert F(64*192,64**2)==3
print("minimum squared 4; Gram determinant 3; center density 1/sqrt(3).")
