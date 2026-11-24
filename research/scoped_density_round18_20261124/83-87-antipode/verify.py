"""Exact arithmetic for the scoped covolume/diameter comparison."""
from fractions import Fraction as F
from itertools import combinations
from math import prod

# Unit-ball volumes in dimensions 0..5, rounded strictly upward using pi<22/7.
p = F(22, 7)
V_up = [F(1), F(2), p, 4*p/3, p*p/2, 8*p*p/15]
# Gram-Schmidt widths <= sqrt((4/3)^i/8), i=0..4.
w = tuple(F(x, 10**7) for x in (3535534,4082483,4714046,5443311,6285394))
assert all(x*x > F(4,3)**j/8 for j,x in enumerate(w))

def coeffs(k):
    return tuple(V_up[j]/2**j*sum(map(prod,combinations(w[:k],k-j))) for j in range(k+1))

def poly(c,x):
    v = F(0)
    for a in c[::-1]:
        v = a+x*v
    return v

def prove(n, upper, uniform_extra=F(0)):
    c = coeffs(88-n)
    # For x=sqrt(beta) in [i/h,(i+1)/h], P(x) <= P((i+1)/h).
    # The n-th-power factor decreases, so evaluate it at i/h; square to
    # remove any square roots even when n is odd.
    h = 1000
    for i in range(3000):
        right = F(i+1,h)
        z = max(F(0), 1+(uniform_extra-F(i,h)**2)/8)
        assert poly(c,right)**2*z**n < upper**2, (n,i)
    print('PASS',n,upper,'uniform source norm surplus allowance',uniform_extra)

for n,x in {83:158,84:169,85:202,86:277,87:450}.items():
    prove(n,F(x,1000))
# Even a uniform extra 1/6 above the parent minimum in all cross-layer pairs
# does not suffice via the same scalar-only argument for ranks 83--85.
for n,x in {83:380,84:417,85:497}.items():
    prove(n,F(x,1000),F(1,6))

# Dimension 83--85 pinned benchmark delta / 2^(n/2) = 1/2.
for n,beta,N,D,ratio in (
    (86,F(1222,2703),7,375708907154061116487005184,F(9,10)),
    (87,F(6417521,14204400),6,115981425005060324321740800,F(19,10)),
):
    # These are the exact center-density data of the cited MW92 packings.
    sq = N*N*((16-beta)/8)**n/D
    assert sq > ratio*ratio
print('PASS exact pinned benchmark comparisons')
