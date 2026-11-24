from fractions import Fraction as F
q=F(376,23)
v=4*(q/18)**48/(q/47)**2
assert q-16==q/47 and q<F(384,23)<18
assert v<F(1,3)
U4=F(636107333,10**9)
# pi>3.14159, e.g. Machin's identity with the alternating series.
def atan_bounds(x,n=20):
    l=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(n))
    return l,l+x**(2*n+1)/(2*n+1)
x,_=atan_bounds(F(1,5))
_,y=atan_bounds(F(1,239))
ell=F(314159,100000)
assert 16*x-4*y>ell
assert 16*U4/ell**2<F(129,125)
print('exact crossover',q,'ratio',float(v))
print('certified full-coset cap relative to D4 <',F(129,125))
from math import factorial
_,U5=atan_bounds(F(1,5))
L239,_=atan_bounds(F(1,239))
pi_hi=16*U5-4*L239
base_occupied_hi=(F(9,2)*pi_hi)**48/(2*47**9*factorial(48))
live_lower=F(6250833919802144,10**37)
assert F(129,125)*base_occupied_hi<live_lower
print('full-coset occupied upper <',float(F(129,125)*base_occupied_hi))
