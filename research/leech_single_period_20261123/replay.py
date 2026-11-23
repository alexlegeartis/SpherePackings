from fractions import Fraction as F
def ratio_sq(k,u):
    return (u/8)**(24+k)*(4/(u-4))**k
for k in range(1,13):
    assert ratio_sq(k,F(5))<1
    assert ratio_sq(k,F(8))==1
for k in range(1,4):
    assert ratio_sq(k,F(257,64))<1
assert ratio_sq(8,F(35,8))<1
assert ratio_sq(9,F(9,2))<1
# Squared projected Leech bound factor before the twelfth power.
def R48(u):return u*u/(18*(u-4))
assert R48(F(6))==1 and R48(F(8))==F(8,9)
# The Hermite constant to the k-th power, then square of the center density.
h=[F(1),F(4,3),F(2),F(4),F(8),F(64,3),F(64),F(256),F(512)]
assert [x/2**i for i,x in enumerate(h,1)]==[F(1,2),F(1,3),F(1,4),F(1,4),F(1,4),F(1,3),F(1,2),F(1),F(1)]
print('PASS exact endpoint inequalities and center-density comparisons')
