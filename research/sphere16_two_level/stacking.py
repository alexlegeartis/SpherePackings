#!/usr/bin/env python3
"""Exact endpoint certificate for quarter-phase stacks with contraction."""
from fractions import Fraction as F
from math import prod
from inputs import TABLE
from verify import SCALE

def wins(n,p,d,Q):
    T=1-1/Q
    m=n//12+2
    c=[1-F(TABLE[n][n-j-1][0],SCALE) for j in range(1,p+1)]
    lim=min(c[-1],F(24-d,12*m)/T)
    for h in sorted({x for x in c if 0<x<=lim}|{lim}):
        need=(1-F(m,2)*T*h)**d*prod(1-T*min(1,h/x) for x in c)
        if (1-T*h)**(n+d)>=need:return False
    return True


STACKING={
 48: ('1150','45','14','7.37','4.74','3.38','2.55','1.98','1.58','1.27','1.04'),
 72: ('2920','82','22.7','11.24','6.88','4.76','3.53','2.75','2.19','1.79','1.48','1.23'),
}

if __name__=='__main__':
    for n,vmax in ((48,5),(72,4)):
        for p,z in enumerate(STACKING[n],1):
            assert all(wins(n,p,d,F(z)) for d in range(1,vmax+1)),(n,p,z)
        print('PASS quarter-phase contractions',n,'additional dimensions through',vmax,STACKING[n])
