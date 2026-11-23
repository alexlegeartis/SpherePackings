"""Directed fixed-point interval proof, involving integer arithmetic only."""
from fractions import Fraction as F
from collections import Counter

DEN=10**12

def ceildiv(a,b):
    return -((-a)//b)

class I:
    """Raw internal endpoints are scaled by DEN; all arithmetic rounds outward."""
    def __init__(self,x,y=None,raw=False):
        if isinstance(x,I):
            assert y is None
            self.l,self.h=x.l,x.h
        elif raw:
            self.l,self.h=x,y
        else:
            x1,y1=F(x),F(x if y is None else y)
            self.l=(x1.numerator*DEN)//x1.denominator
            self.h=ceildiv(y1.numerator*DEN,y1.denominator)
        assert self.l<=self.h
    def __add__(s,t):
        t=I(t);return I(s.l+t.l,s.h+t.h,True)
    __radd__=__add__
    def __neg__(s):return I(-s.h,-s.l,True)
    def __sub__(s,t):return s+-I(t)
    def __rsub__(s,t):return I(t)+-s
    def __mul__(s,t):
        t=I(t); v=(s.l*t.l,s.l*t.h,s.h*t.l,s.h*t.h)
        return I(min(v)//DEN,ceildiv(max(v),DEN),True)
    __rmul__=__mul__
    def __truediv__(s,t):
        t=I(t); assert not t.l<=0<=t.h
        return s*I(DEN*DEN//t.h,ceildiv(DEN*DEN,t.l),True)
    def sq(s):
        return I(0 if s.l<=0<=s.h else min(s.l*s.l,s.h*s.h)//DEN,ceildiv(max(s.l*s.l,s.h*s.h),DEN),True)
    def split(s):
        m=(s.l+s.h)//2
        assert s.l<m<s.h
        return I(s.l,m,True),I(m,s.h,True)

def verify(Q,n,coeff):
    k=I(Q-1)
    c0,c1,c2,d0,d1,d2=(I(x) for x in coeff)
    den=(n+4)*(n+6)
    c3,c4a,c4b,c5a,c5b=map(I,(F(4,n+4),F(24,den),F(6*n,den),F(192,den*(n+8)),F(80*n,den*(n+8))))
    m2coeff=I(F(2,n+2))
    N=I(F(n+2,2))
    one=I(1)
    half=I(F(1,2)); third=I(F(1,3)); fourth=I(F(1,4)); fifth=I(F(1,5))
    todo=[(I(0,1),I(0,1),0)]
    counts=Counter()
    while todo:
        p,u,depth=todo.pop();counts['checks']+=1
        b=one-p
        w=u*p*b
        L=k/(one+k*p)
        # A=k*p/(1+k*p) is monotone in p,k; form endpoints as point intervals.
        slow=I(p.l,p.l,True); shigh=I(p.h,p.h,True)
        A=I((k*slow/(one+k*slow)).l,(k*shigh/(one+k*shigh)).h,True)
        R=L*(half+third*A+fourth*A.sq()+fifth*A.sq()*A/(one-A))
        r=R*w
        lo=-p+u*b;hi=b-u*p
        if (R*(r+c3*hi)).h < m2coeff.l:
            counts['degree3']+=1
            continue
        c=c0+c1*p+c2*w
        d=-w*(d0+d1*p+d2*w)
        Z4=c5a*b+(r-2*c)*c4a
        Z3=c5b*w+(c.sq()-2*d-2*c*r)*c3
        K=(r-2*c)*c4b*w+2*c*d+r*(c.sq()-2*d)+R*d.sq()*N
        if Z4.h<0:
            # Maximize the concave quadratic Z4.h z^2 + t z on [lo.l,hi.h],
            # for both allowed interval extremes t. Products below use the
            # exact integers represented by fixed-point interval endpoints.
            base=(K+Z4*w).h
            ok=True
            for t in (Z3.l,Z3.h):
                a=Z4.h
                zden=-2*a
                znum=t*DEN
                if znum < lo.l*zden:
                    z=lo.l
                elif znum > hi.h*zden:
                    z=hi.h
                else:
                    # Max at interior vertex. Check base / DEN - (t/DEN)^2/(4*a/DEN).
                    ok &= (-4*a)*base+t*t<0
                    continue
                ok &= base*DEN*DEN+a*z*z+t*z*DEN<0
        else:
            z=max(p.h,b.h)
            upper=K.h*DEN*DEN+Z4.h*z*z+max(g*j for g in (Z3.l,Z3.h) for j in (lo.l,hi.h))*DEN
            ok=upper<0
        if ok:
            counts['degree5']+=1
            continue
        assert depth<40,(n,Q,p.l,p.h,u.l,u.h,dict(counts))
        counts['depth']=max(counts['depth'],depth+1)
        if depth%2==0:
            todo.extend((x,u,depth+1) for x in p.split())
        else:
            todo.extend((p,x,depth+1) for x in u.split())
    assert counts['checks']==2*(counts['degree3']+counts['degree5'])-1
    print('PASS',n,Q,dict(counts),flush=True)

if __name__=='__main__':
    import sys
    P=lambda v:tuple(F(x,10000) for x in v)
    if len(sys.argv)==1 or '48' in sys.argv:
        verify(F(2),48,P((1575,5517,19024,5375,10270,-28895)))
    if len(sys.argv)==1 or '72' in sys.argv:
        verify(F(19,10),72,P((1929,-1924,11708,2775,402,-10479)))
