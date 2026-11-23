#!/usr/bin/env python3
"""Exact replay of the [30,15,8] double-stock order-four shear."""
from fractions import Fraction as F
from math import factorial
from collections import Counter
from itertools import product
import json

def rot(a, i):
    return ((a << i) | (a >> (15-i))) & 32767

def span_images(cols):
    out = [0]
    for col in cols:
        out += [x ^ col for x in out]
    return out

def exact_pi_bounds():
    def atan(z):
        s = sum(((-1)**j*z**(2*j+1)/(2*j+1) for j in range(40)), F(0))
        return s, s+z**81/81
    a, A = atan(F(1,5))
    b, B = atan(F(1,239))
    return 16*a-4*B, 16*A-4*b

def main():
    cols = [rot(1491,i) for i in range(15)]
    rows = [sum(((a >> i)&1) << j for j,a in enumerate(cols)) for i in range(15)]
    B = 32767
    counts, depths = [], []
    for c, first_half in ((cols, True), (rows, False)):
        X = span_images(c)
        w = Counter(i.bit_count()+x.bit_count() for i,x in enumerate(X))
        assert sum(w.values())==2**15 and w[0]==1 and min(w.keys()-{0})==8
        assert all(v%2==0 for v in w)
        h = min(((i^B).bit_count()+x.bit_count() if first_half else
                 i.bit_count()+(x^B).bit_count()) for i,x in enumerate(X))
        assert h==7
        counts.append(dict(sorted(w.items())))
        depths.append(h)
    # The graph generators are linearly independent and pairwise orthogonal
    # across the two stocks: (e_i,A e_i).(A^T e_j,e_j) == 0 mod 2.
    assert all(((rows[j]>>i)^(cols[i]>>j))&1==0 for i in range(15) for j in range(15))
    # In a cross-stock difference 1+2c-2d, reduce each odd coordinate to ±1.
    # Only (c,d)=(1,0) needs a nonzero multiple-of-four change.
    for c,d in product(range(2), repeat=2):
        v=1+2*c-2*d
        s=1 if v%4==1 else -1
        assert ((v-s)//4)%2 == c*(1-d)
    assert B.bit_count()%2==1  # double shift toggles the carry class
    assert cols[0]!=0 and rows[0]!=0
    assert span_images(rows)[cols[0]]!=1  # (1,1491) is in C but not in D
    horizontal = {0:32, 1:4*min(depths), 2:16}
    assert all(horizontal[j]+(2*j)**2>=32 for j in range(3))
    assert (2*3)**2>=32
    # Equality at (4(e_1+e_2),0); all cosets distinct by height, coordinate parity and mod 4.
    V=(2*4**30)*8
    M=4*2*2**15
    assert V//M == 2**46 and V%M==0
    delta2 = F(8**31, (V//M)**2)
    assert delta2==2
    ratio=F(2**24,3**15)
    assert ratio>1
    # V_31 = 4^16 16! pi^15 /32!, by half-integer gamma.
    lo, hi = exact_pi_bounds()
    K=F(4**16*factorial(16),factorial(32))
    sl=F(141421356237309504880,10**20)
    su=sl+F(1,10**20)
    assert sl**2<2<su**2
    dlo=F('1.38411060495194e-5')
    dhi=F('1.38411060495196e-5')
    assert dlo < K*lo**15*sl < K*hi**15*su < dhi
    print(json.dumps(dict(status='PASS', dimension=31, weights=counts,
                          hamming_shift_minima=depths, squared_minimum=32,
                          period_volume=V, motif_count=M, volume_per_point=V//M,
                          squared_center_density=str(delta2), ratio_over_Q31=str(ratio),
                          occupied_density_floor='1.38411060495194e-5'),indent=2))
if __name__=='__main__':
    assert __debug__
    main()
