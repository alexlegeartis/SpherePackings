"""Exact cusp polynomial and tail replay; standard library only."""
import json
from fractions import Fraction as F
from pathlib import Path
p=47
one=(1,0)
def mul(u,v):
    a,b=u;c,d=v
    return ((a*c-b*d)%p,(a*d+b*c)%p)
def add(u,v):
    return tuple((a+b)%p for a,b in zip(u,v))
def neg(u):
    return tuple((-a)%p for a in u)
def pw(x,n):
    r=one
    while n:
        if n&1:r=mul(r,x)
        x=mul(x,x);n//=2
    return r
def div(x,y):
    assert y!=(0,0)
    return mul(x,pw(y,p*p-2))
def fm(A,B):
    out=[(0,0)]*(len(A)+len(B)-1)
    for i,a in enumerate(A):
        for j,b in enumerate(B):
            out[i+j]=add(out[i+j],mul(a,b))
    return out
def fp(q):
    return [(x%p,x//p)for x in q]
def slope(x):
    a,b=x
    assert a or b
    return (b*pow(a,-1,p))%p if a else p
def conj(x):
    a,b=x
    return a,-b%p
t0=(18,2);k=(9,9);i=(0,1);a=(5,0);lam=(37,40)
assert pw(t0,48)==neg(one)
assert pw(k,6)==(0,45)
assert pw(lam,6)==a
dilation=div(lam,pw(a,9));assert dilation==(35,1)
assert pw(dilation,6)==(30,0) and pw(dilation,12)==(7,0)
assert pw(pw(dilation,6),23)==(46,0)
change=mul(k,mul((0,2),t0))
w=json.loads(Path(__file__).with_name('witness_polys.json').read_text())
assert set(w)=={'0','1','2','3','16'}
total=set()
for n,obj in w.items():
    assert obj['t0']=='2*z + 18' and obj['k']=='9*z + 9'
    X=fp(obj['X']);Y=fp(obj['Y'])
    assert len(X)==19 and len(Y)==28 and X[-1]!=(0,0) and Y[-1]!=(0,0)
    C=fm(fm(X,X),X);D=fm(Y,Y)
    for j,(c,d) in enumerate(zip(C,D)):
        assert add(d,neg(c))==((int(j==47)-int(j==1))%p,0)
    r=div(X[-1],Y[-1]);original=div(r,change)
    assert [0,0,*original]==obj['jet'] and slope(original)==int(n)
    orb={slope(div(mul(z,pw(dilation,j)),change)) for j in range(6) for z in (r,conj(r))}
    assert not (total&orb)
    total|=orb
    print('literal line',n,'original derivative',original,'orbit',sorted(orb))
assert total==set(range(48))
# Tail certificate for the candidate being tested.
phase=(48,2,4,10)
def depth(m):
    w=[m*a for a in phase]; b=[(v+23)//47 for v in w];r=[v-47*u for v,u in zip(w,b)]
    return sum(a*a for a in r)+(47**2-94*max(map(abs,r)) if sum(b)%2 else 0)
assert sum(phase)%2==0
vals=[depth(m)for m in range(1,47)]
assert min(vals)==478 and [j+1 for j,k in enumerate(vals)if k==478]==[3,44]
assert F(10*min(vals),47**2)>2
print('PASS: every derivative line has an exact squared-norm-18 witness.')

D=set()
for obj in w.values():
    r=div(fp(obj['X'])[-1],fp(obj['Y'])[-1])
    for j in range(276):
        for z in (r,conj(r)):
            D.add(div(mul(z,pw(dilation,j)),change))
assert len(D)==p*p-1 and (0,0) not in D
assert F(16,47**2)*F(10,9)**48/F(314,100)**2<F(116,1000)
print('PASS: all 2,208 nonzero derivatives; arbitrary coloured tail ratio <0.116.')

h20=json.loads(Path(__file__).with_name('H20polys.json').read_text())['H20']
assert h20['t0']=='2*z + 18' and h20['k']=='9*z + 9'
X,Y=fp(h20['X']),fp(h20['Y'])
assert len(X)==21 and len(Y)==31 and X[-1]!=(0,0) and Y[-1]!=(0,0)
C,D=fm(fm(X,X),X),fm(Y,Y)
assert all(add(d,neg(c))==((int(j==47)-int(j==1))%p,0)for j,(c,d) in enumerate(zip(C,D)))
qs=F(18*p,p-1)
assert F(64,1)*F(qs/18)**48/(F(314,100)**2*qs**2)<F(116,1000)
print('PASS: H has a literal norm-20 vector; both endpoint ratios <0.116.')
