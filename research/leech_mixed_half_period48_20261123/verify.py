from fractions import Fraction as F
t=F(371,60)
s=F(11,2)
assert 4<s<t<8
R=[F(8,9)**12,F(11,12)**12,(t/6)**12/2,(s*s/(24*(s-4)))**12,(t*t/(18*(t-4)))**12]
for i,x in enumerate(R,1):
    assert x < F(18,25)
    print('PASS bound',i,x)
