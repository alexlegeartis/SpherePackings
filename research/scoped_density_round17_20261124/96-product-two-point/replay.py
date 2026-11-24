"""Exact, standard-library replay of a formal all-degree two-point witness."""
from itertools import product
from fractions import Fraction as F
from pathlib import Path
from math import comb,isqrt
import json

# Small Gram computations and an exact length cutoff.
def inverse(X):
    n=len(X)
    X=[[F(a) for a in row]+[F(i==j) for j in range(n)] for i,row in enumerate(X)]
    for i in range(n):
        # The matrix is symmetric and all elimination pivots positive,
        # which also certifies positive definiteness.
        assert X[i][i]>0
        X[i]=[a/X[i][i] for a in X[i]]
        for j in range(n):
            if j!=i:
                X[j]=[a-X[j][i]*b for a,b in zip(X[j],X[i])]
    return [x[n:] for x in X]

def span_ok(p,b):
    assert 0<=p<=8 and 0<=b and (p-b)%2==0
    c,d=(p+b)//2,(3*p+b)//2
    G=[[12,6,p,c],[6,18,c,d],[p,c,12,6],[c,d,6,18]]
    H=inverse(G)
    # Any lattice vector with q<=16 has |coordinate_i|^2<=16(G^-1)_ii.
    N=[isqrt(int(16*H[i][i])) for i in range(4)]
    reached=set()
    for v in product(*(range(-n,n+1) for n in N)):
        if not any(v):continue
        GV=[sum(a*z for a,z in zip(row,v)) for row in G]
        q=sum(a*z for a,z in zip(GV,v))
        assert q>=8
        if q<=16 and (v[2] or v[3]):
            w=[h for i in (0,2) for h in (2*v[i+1]-v[i],2*v[i]+v[i+1])]
            s=abs(sum(a*z for a,z in zip(GV,w)))
            if q<=10:assert s<=5*(q-8)
            reached.add((q,s//2))
    return reached

# Fixed-point intervals with exact outward rounding for both signs.
S=10**100
def iv(x):
    x=F(x)*S
    n,d=x.numerator,x.denominator
    return n//d,-((-n)//d)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[1],-x[0]
def sub(x,y):return add(x,neg(y))
def mul(x,y):
    t=[a*b for a,b in product(x,y)]
    return min(t)//S,-((-max(t))//S)
def divi(x,n):
    assert n>0
    return x[0]//n,-((-x[1])//n)
def sc(x,k):
    return (x[0]*k,x[1]*k) if k>=0 else (x[1]*k,x[0]*k)
ONE=iv(1)
s5lo=isqrt(5*S*S)
s5=(s5lo,s5lo+1)
def angles(p,b):
    return add(iv(F(p,12)),mul(iv(F(b,60)),s5)),sub(iv(F(p,12)),mul(iv(F(b,60)),s5))
def powers(x,n):
    y=[ONE]
    for _ in range(n):y.append(mul(y[-1],x))
    return y
def normed(t,n):
    r=[ONE,t]
    for k in range(2,n):
        r.append(divi(sub(sc(mul(t,r[-1]),48+2*k-4),sc(r[-2],k-1)),48+k-3))
    return r

K=100
def envelope(t):
    # Gegenbauer integral bound at k=K; exact moments of u with
    # density proportional to (1-u^2)^22 on [-1,1].
    v=mul(t,t)
    if v[0]<0:v=(0,v[1])
    a=powers(v,K//2);b=powers(sub(ONE,v),K//2)
    z=iv(0);h=F(1)
    for j in range(K//2+1):
        if j:h*=F(2*j-1,2*j+45)
        z=add(z,mul(iv(comb(K//2,j)*h),mul(a[K//2-j],b[j])))
    return z

# Coefficients of the rank-96 Hilbert series, for each fixed sign.
def size(A,q,s):
    return {
      8:{0:A},
      10:{1:231390432000-20*A,3:49886928000-30*A,
          5:2*A+1655821440},
      12:{0:3212*A+1218201435360000,2:726542605152000-770*A,
          4:1240*A+143169114816000,6:185*A+7191815904000,
          8:49886928000-30*A,10:A},
      14:{1:1331791267165632000-20520*A,3:618651644974272000-35720*A,
          5:4620*A+123081209203680000,7:8636968424448000-16480*A,
          9:1240*A+143169114816000,11:231390432000-20*A},
      16:{0:3224640*A+680642990590923360000,2:509650967828070432000-1839820*A,
          4:1840910*A+209200292309552544000,6:43731241296017184000-509740*A,
          8:341440*A+4017806237491956000,10:4620*A+123081209203680000,
          12:726542605152000-770*A,14:231390432000-20*A}
    }[q].get(s,0)

def main():
    C=json.loads(Path(__file__).with_name('counts.json').read_text())
    LA=(0,1662897600)
    seen={}
    for (label,(N,T)),A in zip(C.items(),LA):
        assert N==1218201435360000+3212*A
        W={(p,b):F(w) for p,b,w in T}
        assert len(W)==len(T)
        assert all(w>0 for w in W.values())
        assert 2+sum(W.values())==N
        assert 288+sum(w*p*p for (p,b),w in W.items())==F(3,2)*N
        assert sum(w*b*b for (p,b),w in W.items())==F(15,2)*N
        assert W.get((8,0),0)<=2*753000
        table=[]
        tail=iv(0)
        for (p,b),w in W.items():
            if (p,b) not in seen:
                seen[p,b]=span_ok(p,b)
            orbit=(2 if p else 1)*(2 if b else 1)
            for q,s in seen[p,b]:assert w/orbit<=size(A,q,s),(label,(p,b),(q,s),w,orbit,size(A,q,s))
            u,v=angles(p,b)
            assert all(max(abs(x) for x in y)<S for y in (u,v))
            wt=iv(w/2)
            table.append((wt,normed(u,K),normed(v,K)))
            tail=add(tail,mul(wt,add(envelope(u),envelope(v))))
        assert tail[1]<2*S
        low=2*S
        for i in range(K):
            for j in range(i+1):
                if (i+j)%2 or (i,j) in ((0,0),(1,1),(2,0)):continue
                x=iv(2)
                for w,r,t in table:
                    x=add(x,mul(w,add(mul(r[i],t[j]),mul(r[j],t[i]))))
                assert x[0]>0,(label,i,j,F(x[0],S),F(x[1],S))
                low=min(low,x[0])
        print('PASS',label,N,'support',len(W),'minimum finite moment',float(F(low,S)),
              'all-degree tail bound',float(F(tail[1],S)),
              'mean local p8 per sign',float(W.get((8,0),0)/2))
if __name__=='__main__':main()
