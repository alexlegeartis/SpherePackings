from fractions import Fraction as F
H=F(11,4000)
C=F(129);b=F(19,5)
Q=max(128*b/(3*C),(8/C+4/b)**2)
assert 2**48*H**4 <16800
assert 2**24*(32/(3*C))**4 <16800
assert 2**48*(Q/64)**12/H**4 <16800
print("PASS exact negative-tail cap at c>=129")
z=(1+2*H/7)/15
assert z<F(1,4)
assert 2**48*(2*z)**24/H**4<16800
assert 2**24/F(77,8)**4<16800
print("PASS exact negative-tail cap at b>=7")
from math import factorial
B=F(42*factorial(20),2**21*3**20)
assert B<16800
print("PASS exact moderate-heavy bound at c>=4 and D>=cF/2")
