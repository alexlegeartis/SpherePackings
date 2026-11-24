"""Exact certificate for all bounded metrics with <=3 distinct eigenvalues."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb
from verify import det3

PREC=10**9

def poch(a):
    out=[F(1)]
    for k in range(5):
        out.append(out[-1]*(a+k))
    return out

@lru_cache(None)
def coeffs(n,r,s):
    u,v,w=[poch(F(x,2)) for x in (r,s,n)]
    return [[comb(k,j)*u[k-j]*v[j]/w[k] for j in range(k+1)] for k in range(6)]

@lru_cache(None)
def limit(n,r,s,Q,t):
    target=Q**r*(1+(Q-1)*t)**s
    a,b=(Q-1).as_integer_ratio()
    S=b*PREC
    rhs=target.numerator*S**n
    i,j=0,PREC
    while j-i>1:
        v=(i+j)//2
        if (S+a*v)**n*target.denominator<rhs:
            i=v
        else:
            j=v
    return F(i,PREC)

def test(n,r,s,Q,lo,hi):
    u=limit(n,r,s,Q,lo)
    p=[hi**k for k in range(6)]
    m=[sum(a*b for a,b in zip(c,p)) for c in coeffs(n,r,s)]
    A=[[m[i+j+1]-u*m[i+j] for j in range(3)] for i in range(3)]
    return any(A[i][i]<0 for i in range(3)) or any(A[i][i]*A[j][j]<A[i][j]**2 for i,j in combinations(range(3),2)) or det3(A)<0

def verify(n,Q):
    num=leaves=dep=0
    for r in range(1,n-1):
        for s in range(1,n-r):
            todo=[(F(0),F(1),0)]
            while todo:
                lo,hi,d=todo.pop()
                num+=1
                dep=max(dep,d)
                if test(n,r,s,Q,lo,hi):
                    leaves+=1
                    continue
                assert d<28,(n,r,s,Q,lo,hi)
                m=(lo+hi)/2
                todo.extend(((lo,m,d+1),(m,hi,d+1)))
            # Frozen certificates are explicit rational checks, not global caches.
            limit.cache_clear()
    assert num==2*leaves-(n-1)*(n-2)//2
    print('PASS all three-level shapes',n,Q,'intervals',num,'leaves',leaves,'maxdepth',dep,flush=True)

if __name__=='__main__':
    import sys
    for n,Q in [(48,F(15,4)),(72,F(3))]:
        if len(sys.argv)==1 or str(n) in sys.argv:
            verify(n,Q)
