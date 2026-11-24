"""Exact D5/3D5 geometry, Hoffman spectrum, and radius/density checks."""
from itertools import product
from math import comb
from fractions import Fraction as F
def parity_depth(a):
    # For zero, minimal depth zero. In nonzero residue, center each entry
    # by +/-1; odd sum requires at least one +/-3 correction.
    w=sum(t!=0 for t in a)
    return F(4*(w+3*(w%2)),9)
H=list(product(range(-1,2),repeat=5))
for a in H:
    moves=(k for k in product(range(-1,2),repeat=5) if sum(x+3*y for x,y in zip(a,k))%2==0)
    assert F(4,9)*min(sum((x+3*y)**2 for x,y in zip(a,k)) for k in moves)==parity_depth(a)
assert sorted(set(parity_depth(a)for a in H))==[0,F(8,9),F(16,9),F(24,9),F(32,9)]
# A frequency of weight j contributes 2 in a zero coordinate and -1 otherwise.
spectrum=[]
for j in range(6):
    m=[2]*(5-j)+[-1]*j
    spectrum.append(sum(m[a]*m[b]for a in range(5)for b in range(a+1,5)))
assert spectrum==[40,16,1,-5,-2,10]
assert F(243*5,45)==27
assert F(62,9)+F(8,9)==F(70,9)
assert F(81,32)**2*F(35,36)**85<1
assert F(243,32)**2*F(17,18)**85<1
assert F(243,32)**2*F(35,36)**85>1
# A hypothetical improving full linear rank-five stack must use all
# Hamming-weight 2 classes at r=1, depth 62/9, hence Q=2*weight mod3.
for M in (F(20,3),F(58,9),F(56,9)):
    assert M+F(8,9)<=F(68,9)
assert F(52,9)+F(16,9)==F(68,9)
assert F(29,4)+F(1,2)==F(31,4)
print('PASS all exact geometry, spectra, fixed-period palette and radius comparisons')
