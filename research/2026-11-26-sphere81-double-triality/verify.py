"""Exact verifier, no third-party packages: python verify.py motif.json."""
import sys,json,itertools,collections,math
from functools import lru_cache
from fractions import Fraction
X=[tuple(2*((q>>(k+1))&1)+(q>>5)+4*(k==0)*(q&1) for k in range(4)) for q in range(64)]
def norm(x):
    r=[(a+2)%4-2 for a in x]
    parity=(sum((a-b)//4 for a,b in zip(x,r))&1)
    n=sum(a*a for a in r)+parity*min(16-8*abs(a) for a in r)
    assert n%4==0
    return n//4
def minus(a,b):return tuple(x-y for x,y in zip(a,b))
@lru_cache(None)
def enc(x):return next(i for i,y in enumerate(X) if norm(minus(x,y))==0)
D=[[norm(minus(x,y)) for y in X] for x in X]
ADD=[[enc(tuple(a+b for a,b in zip(x,y))) for y in X] for x in X]
NEG=[enc(tuple(-v for v in x)) for x in X]
def add(x,y):return (((x^y)>>12)<<12) + 64*ADD[x>>6&63][y>>6&63]+ADD[x&63][y&63]
def neg(x):return (x>>12<<12) +64*NEG[x>>6&63]+NEG[x&63]
def dist(x,y):return ((x^y)>>12).bit_count()+D[x>>6&63][y>>6&63]+D[x&63][y&63]
def color(q):return 2*(q>>5)+((q>>1&15).bit_count()%2)
def check(A):
    assert len(A)==len(set(A))==288 and all(0<=x<32768 for x in A)
    S=set(A);K=collections.Counter(dist(x,y) for i,x in enumerate(A) for y in A[:i])
    assert min(K)==4 and K[4]==16020
    def common(t):return sum(add(x,t) in S for x in A)
    H=sorted(t for y in A for t in [add(y,neg(A[0]))] if common(t)==len(A))
    assert H==[0]
    print('n, contact pairs, stabilizer:',len(A),K[4],H)
    cols={(color(x>>6&63),color(x&63)) for x in A};print('affine colors',sorted(cols))
    d=collections.defaultdict(list)
    for v in itertools.product(range(-4,5),repeat=4):
        n=sum(k*k for k in v)
        if n>16 or n%4 or any((z-v[0])%2 for z in v):continue
        d[n//4].append((v,enc(v)))
    by=collections.defaultdict(list);ov={}
    for o in itertools.product(range(-2,3),repeat=3):
        m=4-sum(i*i for i in o)
        if m<0:continue
        p=sum((v%2)<<j for j,v in enumerate(o))
        for k in range(m+1):
            for a,q in d[k]:
                for b,r in d[m-k]:
                    t=4096*p+64*q+r
                    if t not in ov:ov[t]=common(t)
                    by[ov[t]].append(tuple(2*v for v in o)+a+b)
    piv={}; cut={}
    for v in sorted(by,reverse=True):
        for x in by[v]:
            z=list(x)
            for i in sorted(piv):
                w=piv[i]
                if z[i]:
                    a,b=w[i],z[i]; z=[a*k-b*j for k,j in zip(z,w)]
                    g=math.gcd(*z)
                    if g>1:z=[k//g for k in z]
            for i,k in enumerate(z):
                if k:
                    piv[i]=z;break
        for r in (7,11):
            if r not in cut and len(piv)>=r:cut[r]=v
        if len(piv)==11:break
    print('overlap rate ceilings for seven and eleven independent length-two periods:',cut)
    assert len(piv)==11 and cut=={7:272,11:272}
    # The above ceilings are upper bounds, not necessarily attained: any
    # subpacking with given periods has intensity bounded by its overlap
    # with its translate by EVERY one of these periods. A translation
    # outside Z^3 x D4* x D4* has zero overlap.
if __name__=='__main__':check(json.load(open(sys.argv[1]))['motif'])
def charts():
    H=[0,1,30,31];assert {ADD[x][y] for x in H for y in H}==set(H)
    assert {tuple(sorted(D[q][h] for h in H)) for q in range(64)}=={(0,4,4,4),(1,1,3,3),(2,2,2,2)}
    A=[p*4096+64*x+ADD[x][h] for p,h in zip([0,3,5,6],H) for x in H]
    assert min(2*((x^y)>>12).bit_count()+D[x>>6&63][y>>6&63]+D[x&63][y&63] for i,x in enumerate(A) for y in A[:i])>=8
    def omega(a,b):return ((a&(((b&5)<<1)|((b&10)>>1))).bit_count()&1)
    planes={tuple(sorted((0,u,v,u^v))) for u in range(1,16) for v in range(u+1,16) if not omega(u,v)}
    assert len(planes)==15 and sum(any((q>>2)==0 and q for q in p) for p in planes)==9
    for q in range(64):
        for r in range(64):
            assert (sum(x*y for x,y in zip(X[q],X[r]))%4)//2 == omega(color(q),color(r))
    print('All endpoint charts and the 16-site witness, 15 symplectic plane types and scalar-product pairing verified.')
if __name__=='__main__':charts()
