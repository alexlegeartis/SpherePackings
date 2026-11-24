"""Complete affine golden unit identifications for norm 10,12,14."""
from joint3 import jointrows3
from affine import advance,shift
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
def itinerary(y,m):
    q,p,c,s=y
    rho=q-(p*p+c*c)/10
    if rho==0:
        # In the anchor plane 10 ||y||^2=p^2+c^2. 120 and 140
        # are not sums of squares; at norm 10 only the given four.
        assert (q,p,c,s) in [(10,10,0,-10),(10,-10,0,-10),(10,0,10,10),(10,0,-10,10)]
        return
    assert rho>0
    yield y
    for k in (-1,1):
        z=y
        for _ in range(50):
            z=advance(z,k)
            q,p,c,s=z
            h=q-(p*p+c*c)/10
            t=s-(-p*p+4*p*c+c*c)/10
            if h>m and k*t>=0:break
            yield z
        else:raise ArithmeticError(y)
def relations(parts,qn=(10,12,14)):
    Z=[(d,bb) for d,B in zip(qn,parts) for bb in B]
    lut={(d,*y):i for i,(d,bb) in enumerate(Z) for y in bb}
    assert len(lut)==sum(len(bb) for d,bb in Z)
    E=set()
    for i,(d,bb) in enumerate(Z):
        for p,c,s in bb:
            y=tuple(map(F,(d,p,c,s)))
            for z in itinerary(y,max(qn)):
                p0,c0=map(int,z[1:3])
                assert z[1:3]==(p0,c0)
                for u,v in product(range((-10-p0+9)//10,(10-p0)//10+1),
                                   range((-10-c0+9)//10,(10-c0)//10+1)):
                    w=shift(z,u,v)
                    if w[0] not in qn:continue
                    assert all(t.denominator==1 for t in w)
                    j=lut.get(tuple(map(int,w)),-1)
                    if i!=j:E.add(tuple(sorted((i,j))))
    return sorted(E)
def classes(m,E):
    P=list(range(m+1))
    def find(i):
        while P[i]!=i:
            P[i]=P[P[i]]
            i=P[i]
        return i
    for i,j in E:
        a,b=sorted((find(i+1),find(j+1)))
        P[b]=a
    return [find(i+1)-1 for i in range(m)]
def load():
    R,V,M,A,B,C=jointrows3()
    E=relations((A,B,C));m=len(A)+len(B)+len(C)
    assert m==1391
    J=classes(m,E)
    return R,V,M,(A,B,C),E,J
def old():
    return {int(k):F(v) for k,v in json.loads(Path(__file__).with_name('joint3-counts.json').read_text()).items()}
if __name__=='__main__':
    R,V,M,P,E,J=load()
    T=old()
    fail=[(i,j,T.get(i,0),T.get(j,0)) for i,j in E if T.get(i,0)!=T.get(j,0)]
    print('relations',len(E),'classes',len(set(J)-{-1}),'forced zero',J.count(-1),flush=True)
    print('Failures',len(fail),fail[:12])
    print('First failure labels',P[0][fail[0][0]] if fail and 0<=fail[0][0]<len(P[0]) else None)
