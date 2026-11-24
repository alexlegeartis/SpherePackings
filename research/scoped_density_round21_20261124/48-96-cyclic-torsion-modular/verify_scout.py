import gzip,json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

def mm(a,b):
    return [[sum(x*y for x,y in zip(r,c))for c in zip(*b)]for r in a]
def tr(a):return list(map(list,zip(*a)))
def det(a):
    a=[list(map(int,r))for r in a];s=1;d=1
    for k in range(len(a)-1):
        i=next((i for i in range(k,len(a)) if a[i][k]),None)
        if i is None:return 0
        if i!=k:s=-s;a[i],a[k]=a[k],a[i]
        b=a[k][k]
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                v=a[i][j]*b-a[i][k]*a[k][j]
                assert v%d==0
                a[i][j]=v//d
            a[i][k]=0
        d=b
    return s*a[-1][-1]
J=json.load(gzip.open(Path(__file__).with_name("p48_scout.json.gz"),'rt'))
B=J["basis"];P=J["parent"];G=J["G"];assert abs(det(B))==1 and det(P)==1
assert mm(mm(B,P),tr(B))==G
def quadratic(d):return sum(x*sum(y*a for y,a in zip(d,gi))for x,gi in zip(d,G))
counts=Counter()
far={}
for p,g,c,z,q in J["rows"]:
    assert p in (5,6,7,13) and 0<=g<13 and len(c)==48 and all(0<=u<p for u in c)
    if z is None:
        assert q is None and p==6 and all(u in (0,3)for u in c)
        # Half-period bound from the separate half-torus theorem.
    else:
        assert len(z)==48 and all(isinstance(a,int) for a in z)
        assert q==quadratic([a-p*b for a,b in zip(c,z)]) and 0<=q<=4*p*p and q%2==0
        if not (p==6 and all(u in (0,3) for u in c)):assert q<4*p*p
        far[p]=max(far.get(p,F(0)),F(q,p*p))
    counts[p,z is None]+=1
assert counts=={(5,False):832,(6,True):90,(6,False):742,(7,False):832,(13,False):832}
print("Replay actual source transform and all finite scout blockers:",dict(counts))
print("Largest displayed rational squared blocker:",far)
