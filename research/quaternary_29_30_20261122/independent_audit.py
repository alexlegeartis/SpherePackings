from collections import Counter
from fractions import Fraction as F
from itertools import product
from functools import reduce
from operator import xor
N=30
T=tuple(range(2,32))
V=(1620470160,2896865964,3224620080,2695517520,1617532560,4293980400)
W=(2155872256,0,142606336,8421376,545259520,2147516416)
column=lambda j: 32+j
bit=lambda x,j: x>>j&1
maskT=sum(1<<j for j in T)
assert all((v|w)&~maskT==0 for v,w in zip(V,W))
A=tuple(tuple(bit(v,j)+2*bit(w,j) for j in T) for v,w in zip(V,W))
# Enumerate the image of the parity-check map; T_0 has 2^(N-6) words.
image={0}
for j in T:image|={x^column(j) for x in image}
assert len(image)==64
Hamming=lambda seq:reduce(xor,(column(j) for j,x in zip(T,seq) if x),0)
assert all(Hamming([x%2 for x in a])==0 for a in A) # doubles of generators belong to 2T_0
# No weight 1,2,3 word in T_0: first and third odd sizes fail parity,
# while two different columns give unequal labels. Exhibit a weight-four word.
assert len({column(j) for j in T})==N and reduce(xor,map(column,(2,3,4,5)))==0
tables=Counter(); seen=set(); status={}
for f in product((0,1),repeat=6):
 r=tuple(sum(b*a[j] for a,b in zip(A,f))%4 for j in range(N))
 s=tuple(x%2 for x in r)
 assert s not in seen;seen.add(s)
 if not any(f): assert not any(s);continue
 # trellis: state is the Hamming parity check of the added binary vector.
 # Each node considers both t_j choices in the quadrature cost exactly.
 costs=[0]+[99]*63
 for j,a in zip(T,r):
  weight=(0,1,4,1)
  u=column(j);v=weight[a];w=weight[(a+2)%4]
  costs=[min(costs[k]+v,costs[k^u]+w) for k in range(64)]
 d=costs[0]
 assert d>=16,(f,r,d)
 tables[d]+=1;status[sum(b<<i for i,b in enumerate(f))]=d
assert len(seen)==64
card=2**6*2**24
assert F(2)**N/card==1
pair=next((i,j,sum(x*y for x,y in zip(A[i],A[j]))%4) for i in range(6) for j in range(i,6) if sum(x*y for x,y in zip(A[i],A[j]))%4)
print('independent quaternary trellis distribution',dict(sorted(tables.items())), 'covolume 1, exact minimum squared 4')
print('nonintegral pair: modulus-four dot',pair)
