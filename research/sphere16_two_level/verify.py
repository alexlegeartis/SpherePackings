#!/usr/bin/env python3
"""Integer/Fraction-only proof replay for two-level and low-rank metric bounds."""
from fractions import Fraction as F
from itertools import permutations
from math import prod

from inputs import TABLE, CONTRACTION

SCALE=10**12
Q_DEN=1000

def beta_moments(n,r):
    ans=[F(1)]
    for k in range(5):
        ans.append(ans[-1]*F(r+2*k,n+2*k))
    return ans

def det3(a):
    return sum((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
               *prod(a[i][p[i]] for i in range(3)) for p in permutations(range(3)))

def verify():
    assert set(TABLE)==set(CONTRACTION)=={48,72}
    for n,table in TABLE.items():
        assert len(table)==n-1
        for r,(u,q) in enumerate(table,1):
            m=beta_moments(n,r)
            x=F(u,SCALE)
            assert 0<x<F(r,n)
            assert det3([[m[i+j+1]-x*m[i+j] for j in range(3)] for i in range(3)])<0
            Q=F(q,Q_DEN)
            assert Q>1 and (1+(Q-1)*x)**n < Q**r, (n,r,q)
        Q=[F(q,Q_DEN) for _,q in table]
        assert min(Q)>=({48:F(15,4),72:F(3)}[n])
        assert Q[0]>=6 and Q[-1]>=({48:30000,72:38000}[n])
        c=[1-F(table[n-j-1][0],SCALE) for j in range(1,13)]
        assert all(F(0)<a<b<F(1) for a,b in zip(c,c[1:]))
        for k,q in enumerate(CONTRACTION[n],1):
            Q=F(q,Q_DEN)
            assert Q>1
            T=1-1/Q
            for h in c[:k]:
                assert (1-T*h)**n < prod(1-T*min(1,h/ci) for ci in c[:k]),(n,k,q,h)
        print('rank',n,'all two-level ratios up to at least',min(F(q,Q_DEN) for _,q in table))
        print('rank',n,'arbitrary rank-k contraction ratios, k=1..12:',
              ', '.join(str(F(q,Q_DEN)) for q in CONTRACTION[n]))
    print('PASS: exact Beta moments, matrix signs, and rational density comparisons')

if __name__=='__main__':
    verify()
