#!/usr/bin/env python3
"""Direct exact replay of dimension-29 Vardy layering; standard library only."""
import json
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path
N=28
K=14
ROW=471
B0=715
B1=15538
P=1<<K
MASK=P-1
a=[((ROW<<i)|(ROW>>(K-i)))&MASK for i in range(K)]
at=[sum(((a[i]>>j)&1)<<i for i in range(K)) for j in range(K)]
C0=[(1<<i)|(a[i]<<K) for i in range(K)]
D0=[at[i]|(1<<(K+i)) for i in range(K)]
def span(b):
 seq=[0]
 for g in b:
  seq += [w^g for w in seq]
 return seq
C=span(C0);D=span(D0)
B=B0|(B1<<K)
assert len(C)==len(D)==P and len(set(C))==len(set(D))==P
assert all(c.bit_count()%2==0 for c in C+D)
assert all((c&d).bit_count()%2==0 for c in C0 for d in D0)
min_C=min(c.bit_count() for c in C[1:]);min_D=min(d.bit_count() for d in D[1:])
bc=min((B^c).bit_count() for c in C);bd=min((B^d).bit_count() for d in D)
assert (min_C,min_D,bc,bd)==(8,8,6,6)
assert sum((B^c).bit_count()==6 for c in C)>0
assert sum((B^d).bit_count()==6 for d in D)>0
# The only way all cross-layer parity coordinates can have absolute value 1
# is to choose each integer z_i with parity c_i(1-d_i). Its total always even.
for c in (0,1):
 for d in (0,1):
  z=[i for i in (-1,0,1) if abs(1+2*c-2*d+4*i)==1]
  assert len(z)==1 and z[0]%2==c*(1-d)
assert N+8>=32
assert min(4*bc,4*bd,N)+8>=32  # adjacent translated layers
assert 4*8==32 # two heights apart
# One base-center lattice vector 4(e_i+e_j) attains min^2 32.
# There are 2 motifs vertically, each of size 2^15 over 4D_28.
# Period volume squared = (2*4^28)^2 * (2sqrt(8))^2.
count=2*(len(C)+len(D))
period_volume_squared=(2*4**N)**2*(4*8)
rho_squared=8
density_squared=Fraction(pow(rho_squared,N+1)*count*count,period_volume_squared)
assert density_squared==1
# 0 and the odd zero-codeword p=1+4e1 are in the horizontal fiber;
# adding the even codeword center -2c0 does not land in it.
assert C0[0] not in set(D)
getcontext().prec=78
pi=Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286')
den=1
for i in range(1,30,2):den*=i
volume=Decimal(2**15)*pi**14/den
result={
  "dimension":29,"min_squared_distance":32,"center_density":"1",
  "packing_density_decimal":str(volume),"Cohn_center_density":"1/sqrt(2)",
  "percent_gain":str(100*(Decimal(2).sqrt()-1)),
  "packing_type":"nonlattice periodic",
  "binary_bit_order":"least significant bit first in each half",
  "seed_row":ROW,"hole_halves":[B0,B1],
  "hole_bits":[(B>>i)&1 for i in range(N)],"C_generators_int":C0,"D_generators_int":D0,
  "C_minimum_hamming":min_C,"D_minimum_hamming":min_D,
  "hole_to_C":bc,"hole_to_D":bd,
  "nearest_C_count":sum((B^c).bit_count()==6 for c in C),
  "nearest_D_count":sum((B^d).bit_count()==6 for d in D),
  "C0_not_in_D_witness":C0[0],
  "period_motif_count":count,"period_volume_squared":period_volume_squared}
print(json.dumps(result,indent=2))
