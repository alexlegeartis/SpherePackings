"""Exact boxes for arbitrary low-rank positive-semidefinite stretches."""
from fractions import Fraction as F
from itertools import combinations
from functools import lru_cache
from math import factorial, prod
from verify import det3

D=10**9
@lru_cache(None)
def rndmom(n, hi):
    """Moments of sum_i hi[i] X_i^2 for X uniform on the unit sphere."""
    aa=[F(1)]
    nn=[F(1)]
    for k in range(1,6):
        aa.append(aa[-1]*F(2*k-1,2*k))
        nn.append(nn[-1]*F(n+2*k-2,2))
    z=[F(1),0,0,0,0,0]
    for b in (F(1),)+hi:
        c=[v*b**i for i,v in enumerate(aa)]
        z=[sum(z[j]*c[k-j] for j in range(k+1)) for k in range(6)]
    return [z[k]*factorial(k)/nn[k] for k in range(6)]

@lru_cache(None)
def threshold(n,s,lo):
    """Largest rational v=k/D with (1+s*v)^n < prod_i(1+s*lo_i)."""
    target=(1+s)*prod(1+s*t for t in lo)
    num,den=s.numerator,s.denominator
    scale=D*den
    base=target.numerator*scale**n
    j,k=0,D
    assert target < (1+s)**n
    while k-j>1:
        i=(j+k)//2
        if (scale+num*i)**n*target.denominator<base:
            j=i
        else:
            k=i
    return F(j,D)

def cell(n,s,lo,hi):
    x=threshold(n,s,lo)
    m=rndmom(n,hi)
    a=[[m[i+j+1]-x*m[i+j] for j in range(3)]for i in range(3)]
    return any(a[i][i]<0 for i in range(3)) or any(a[i][i]*a[j][j]<a[i][j]**2 for i,j in combinations(range(3),2)) or det3(a)<0

def verify(n,k,Q):
    print('starting stretch',n,k,Q,flush=True)
    todo=[((F(0),)*(k-1),(F(1),)*(k-1),0)]
    num=closed=empty=maxdepth=0
    while todo:
        lo,hi,d=todo.pop()
        num+=1
        maxdepth=max(maxdepth,d)
        lo=tuple(max(lo[i:]) for i in range(k-1))
        hi=tuple(min(hi[:i+1]) for i in range(k-1))
        if any(a>b for a,b in zip(lo,hi)):
            empty+=1
            continue
        if cell(n,Q-1,lo,hi):
            closed+=1
            continue
        assert d<42, (n,k,Q,lo,hi)
        j=max(range(k-1),key=lambda i:hi[i]-lo[i])
        m=(lo[j]+hi[j])/2
        todo.extend(((lo,hi[:j]+(m,)+hi[j+1:],d+1),(lo[:j]+(m,)+lo[j+1:],hi,d+1)))
    assert num==2*(closed+empty)-1
    print('PASS stretch',n,k,Q,'boxes',num,'closed',closed,'empty',empty,'depth',maxdepth,flush=True)
    rndmom.cache_clear()
    threshold.cache_clear()

if __name__=='__main__':
    import sys
    KS={48:{2:F(11,2),3:F(5),4:F(47,10),5:F(9,2),6:F(43,10),7:F(41,10)},72:{2:F(26,5),3:F(19,4),4:F(22,5),5:F(41,10),6:F(39,10),7:F(15,4)}}
    for n,ks in KS.items():
        if len(sys.argv)>1 and str(n) not in sys.argv:
            continue
        for k,q in ks.items():
            verify(n,k,q)
