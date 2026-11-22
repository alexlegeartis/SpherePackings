#!/usr/bin/env python3
"""Exact finite arithmetic for the MW92 cusp–D4 lattice.
The underlying Mordell–Weil rank/covolume and geometric specialization
and height arguments are stated in README.md.
"""
from fractions import Fraction as F
from math import factorial
from itertools import product
from collections import Counter
from decimal import Decimal,localcontext
import json,sys
p=47
def mul(a,b):
 x,y=a;u,v=b
 return ((x*u-y*v)%p,(x*v+y*u-y*v)%p)
def poww(a,n):
 z=(1,0)
 for _ in range(n):z=mul(a,z)
 return z
assert all((x*x+x+1)%p for x in range(p))
t=(3,23); w=(0,1); c=poww(t,24)
assert c==(4,8) and poww(t,48)==(p-1,0)
assert mul(w,c)==((-c[1])%p,(c[0]-c[1])%p)
assert (c[0]*mul(w,c)[1]-c[1]*mul(w,c)[0])%p
assert sum(poww(u,48)==(p-1,0) for u in product(range(p),repeat=2))==48
C=((1,0),(26,33),(36,35),(35,37)) # columns relative to (c,wc)
assert (C[0][0]*C[1][1]-C[0][1]*C[1][0])%p
def phi(coeff):
 return tuple(sum(q*C[j][i] for j,q in enumerate(coeff))%p for i in range(2))
assert phi((1,0,0,0))==(1,0)
count=Counter()
for z in product(range(-4,5),repeat=4):
 n=sum(x*x for x in z)
 if not 0<n<=16 or sum(z)%2:continue
 # D4 simple basis e0-e1,e1-e2,e2-e3,e2+e3.
 a=z[0];b=z[0]+z[1];d=sum(z)//2;c0=d-z[3]
 assert z==(a,b-a,c0+d-b,d-c0)
 assert phi((a,b,c0,d))!=(0,0)
 count[n]+=1
assert count=={2:24,4:24,6:96,8:24,10:144,12:96,14:192,16:24}
assert sum(count.values())==624
size=p**2;covol=2*size*p**7
center=F(9,2)**48/covol
D=F(19889279550909941977,75000000000000000)
current_lattice_center_squared=F(4)**96/(p**14*D)
r=center**2/current_lattice_center_squared
assert r>F(10516,10000)**2
bcs=F(3)**48/2**24
assert center/bcs>F(212417,100000)
def atan_interval(x,m=42):
 q=sum(F((-1)**j,2*j+1)*x**(2*j+1) for j in range(m))
 rem=x**(2*m+1)/(2*m+1)
 return (q,q+rem) if m%2==0 else (q-rem,q)
l,h=atan_interval(F(1,5));L,H=atan_interval(F(1,239))
pi_lo=16*l-4*H;pi_hi=16*h-4*L
lo,hi=[center*x**48/factorial(48) for x in (pi_lo,pi_hi)]
def dec(f):
 with localcontext() as ct:
  ct.prec=30
  return Decimal(f.numerator)/Decimal(f.denominator)
def floor_digits(x,scale):
 y=x*10**scale
 return F(y.numerator//y.denominator,10**scale)
out={'rank':96,'minimum_squared':18,'cusp':[3,23],
     'root_label':[4,8], 'D4_specialization_columns':C,'D4_nonzero_short_counts':dict(sorted(count.items())),
     'covolume':covol,'covolume_formula':'2*47^9','center_density':str(center),
     'center_density_approx':str(dec(center)),
     'occupied_density_floor':str(dec(floor_digits(lo,38))),
     'gain_vs_previous_lattice_percent_approx':str((dec(r).sqrt()-1)*100),
     'gain_vs_BCS96_percent_approx':str((dec(center/bcs)-1)*100),
     'squared_ratio_vs_previous_lattice':str(r),
     'effective_det_at_min16':str(F(4418)**2*F(8,9)**96)}
assert hi-lo<F(1,10**75)
assert hi<F(62508,10**26)   # 6.2508e-22, below known unrestricted
print(json.dumps(out,indent=2))
