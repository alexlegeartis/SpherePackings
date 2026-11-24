"""Exact scalar-series and distance gates. See README for the mathematical proof."""
from fractions import Fraction as Q
from math import ceil,floor
from sympy import Matrix,Rational

def ring(N):
    def const(x):return (x,)+(0,)*(N-1)
    def add(a,b):return tuple(x+y for x,y in zip(a,b))
    def mul(a,b):return tuple(sum(a[j]*b[i-j] for j in range(i+1)) for i in range(N))
    def power(a,n):
        assert n>=0
        out=const(1)
        for _ in range(n):out=mul(out,a)
        return out
    return const,add,mul,power

def anisotropic_count(n):
    c,add,mul,powr=ring(7)
    th3=(1,2,0,0,2,0,0);th4=(1,-2,0,0,2,0,0)
    X=(0,16,0,64,0,96,0);Y=powr(th4,4);W=add(X,Y);Z=add(X,tuple(2*b for b in Y))
    e=n//12-1;d=n//4-1-2*e
    M=[mul(mul(powr(th3,2),powr(th4,2)),mul(powr(Y,e),mul(powr(W,e),mul(powr(Z,d-2*u),powr(X,2*u))))) for u in range(d//2+1)]
    l=n//24
    A=Matrix([[v[2*i] for v in M] for i in range(l+1)])
    z,par=A.gauss_jordan_solve(Matrix([1]+[0]*l))
    assert not len(par)
    return int(4*16**e*sum(z))

Hmax={n:anisotropic_count(n) for n in (48,72)}
assert Hmax=={48:108,72:6272}

def series_case(n,r,H):
    D=n//12;l=n//24
    c,add,mul,powr=ring(l+1)
    th3=(1,2,0,0)[:l+1];th4=(1,-2,0,0)[:l+1]
    X=(0,16,0,64)[:l+1];Y=powr(th4,4);W=add(X,Y)
    P=mul(powr(th4,r),powr(th3,4-r)) if r else c(1)
    i=2*D+(r==0)
    alpha=Q(4-r,4) if r else 0
    j=ceil(H/2-alpha);d=n//4-(r>0)-i-j
    if d<0:return None
    M=[mul(P,mul(powr(Y,i),mul(powr(W,j),mul(powr(X,u),powr(Y,d-u))))) for u in range(d+1)]
    A=Matrix([[v[t] for v in M] for t in range(l+1)])
    try: z,par=A.gauss_jordan_solve(Matrix([1]+[0]*l))
    except ValueError:return None
    c0=Rational(2)**(r+4*i-n//2)*z[-1]
    c1=abs((2**((4-r if r else 0)+4*j))*sum((-1)**(d-u)*z[u] for u in range(d+1)))
    return c0,c1,2*(alpha+j),Q(i,2)+Q(r,8),len(par)

def geom(n,r,H):
    D=n//12
    parity=Q((r-int(2*H))%4,4)
    assert parity in (0,Q(1,2))
    m=parity
    while m*m<H/3:m+=1
    q=max(D-Q(1,2)+H/4-parity/2,D+H/4-m/2)
    return Q(r,8)+Q(floor(2*(q-Q(r,8))),2)

expected={
(48,1,Q(7,2)):(3728,516,Q(7,2),Q(33,8),0),
(48,2,Q(3)):(8448,304,Q(3),Q(17,4),0),
(72,0,Q(4)):(2514944,5840,Q(4),Q(13,2),0),
(72,1,Q(11,2)):(214592,29648,Q(11,2),Q(49,8),0),
(72,3,Q(7,2)):(1112832,10108,Q(9,2),Q(51,8),0),
(72,3,Q(9,2)):(1112832,10108,Q(9,2),Q(51,8),0),
}
for key,entry in expected.items():assert series_case(*key)==entry
shelf={
 (48,1,Q(7,2)):Hmax[48],
 (48,2,Q(3)):96,
 (72,0,Q(4)):144,
 (72,1,Q(11,2)):Hmax[72],
 (72,3,Q(7,2)):1280,
 (72,3,Q(9,2)):1280,
}
a=Q(1,9)
assert 2*72*(1-a*a)/(1-72*a*a)==1280
for key,cap in shelf.items():
    assert expected[key][1]>cap
    print('SIGNED CUSP IMPOSSIBLE',key,'needs',expected[key][1],'available',cap)
remaining={}
for n in (48,72):
    D=n//12
    remain=[]
    for r in range(4):
        for h2 in range(n//24+1,2*D+1):
            H=Q(h2,2)
            if (H-Q(r,2)).denominator>1:continue
            u=geom(n,r,H)
            if u<=D:continue
            key=(n,r,H)
            dat=series_case(*key)
            if dat is None or key in shelf:continue
            lead,_,_,q,p=dat
            if not p and lead>0:u=min(u,q)
            assert u>D
            remain.append((r,H,u))
    remaining[n]=remain
assert remaining=={
48:[(1,Q(5,2),Q(33,8))],
72:[(1,Q(5,2),Q(49,8)),(1,Q(7,2),Q(49,8)),(1,Q(9,2),Q(49,8)),(2,Q(3),Q(25,4)),(2,Q(4),Q(25,4)),(2,Q(5),Q(25,4))]}
assert series_case(48,1,Q(5,2))==(3728,516,Q(7,2),Q(33,8),0)
for n,v in remaining.items():
    print('PRIOR (before new shell gates)',n,'remaining possible above baseline',v)
