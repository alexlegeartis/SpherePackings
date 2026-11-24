"""Rational Delsarte certificate for the closest layer inside q=12,s=0."""
from fractions import Fraction as F
D=94
C={3:F(4671701,12500),5:F(654129257,100000),7:F(532077891,20000),
   9:F(3550223021,2000),10:F(822315292389,100000),11:F(630060429287,50000),
   12:F(10181947194587,100000),13:F(4698872419817,25000)}
def mono(a):
 t=[F(1),a]
 for k in range(2,14):
  t.append(((D+2*k-4)*a*t[-1]-(k-1)*t[-2])/(D+k-3))
 return 1+sum(v*t[i] for i,v in C.items())
assert all(v>0 for v in C.values())
assert all(mono(F(3*j-4,20))<0 for j in range(-4,5))
assert mono(F(1))<313_000_000
print("PASS: q12-to-q8 joint orthogonal local cap below 313 million")
