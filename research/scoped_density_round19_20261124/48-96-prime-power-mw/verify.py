"""Exact finite part of the tame crystalline regulator and short section proof."""
from collections import Counter

def orbits(A,p,M):
    A=set(A)
    while A:
        u=min(A); v=u; S=[]
        while v in A:
            A.remove(v)
            S.append(v)
            v=(p*v)%M
        assert v==u
        yield S

def exponents(p,q,k):
    M=k*(q-1)
    hol={k*i-j for j in range(1,k) for i in range(1,(j*q-1)//k+1)}
    A={a for a in range(M) if a%k}
    assert len(hol)*2==len(A) and set((-a)%M for a in hol)==A-hol
    slope_shapes=Counter()
    for orbit in orbits(A,p,M):
        count=sum(x in hol for x in orbit)
        assert count*2==len(orbit)
        slope_shapes[len(orbit)]+=1
    totals=Counter()
    active={a for a in A if a%k in (1,k-1)}
    for orbit in orbits(active,p,M):
        v=[0]
        for a in orbit:
            v.append(v[-1]+int(a in hol)+int(a%k==k-1)-1)
        assert v[-1]==0
        totals[sum(x-min(v) for x in v[:-1])]+=1
    return len(active),slope_shapes,totals

stats=[]
for p,q,k in [(5,5,6),(5,25,6),(7,7,4),(7,49,4)]:
    n,S,T=exponents(p,q,k)
    sigma=sum(a*m for a,m in T.items())
    print(p,q,k,'rank',n,'H1 orbit counts',dict(S),'sigma summands',dict(T),'sigma',sigma)
    stats.append((n,sigma,T))
assert stats[0]==(8,0,{0:4})
assert stats[1]==(48,24,{3:1,5:1,7:1,9:1})
assert stats[2]==(12,1,{0:5,1:1})
assert stats[3]==(96,48,{1:2,3:4,5:4,7:2})

# Field F7[L]/(L²+2), and the literal quadratic-square test
assert 5 not in {x*x%7 for x in range(7)}
def add(a,b):return ((a[0]+b[0])%7,(a[1]+b[1])%7)
def mul(a,b):return ((a[0]*b[0]+5*a[1]*b[1])%7,(a[0]*b[1]+a[1]*b[0])%7)
def neg(a):return((-a[0])%7,(-a[1])%7)
zero=(0,0); one=(1,0); lam=(0,1)
def pa(A,B):return [add(x,y) for x,y in zip(A+[zero]*(max(len(A),len(B))-len(A)),B+[zero]*(max(len(A),len(B))-len(B)))]
def pm(A,B):
    R=[zero]*(len(A)+len(B)-1)
    for i,a in enumerate(A):
        for j,b in enumerate(B):R[i+j]=add(R[i+j],mul(a,b))
    return R
def poly(roots):
    A=[one]
    for a in roots:A=pm(A,[(-a%7,0),one])
    return A
A=poly(range(4));B=poly(range(4,7))
assert pm(A,B)==[zero,(6,0)]+[zero]*5+[one]
C=pa(A,[mul(lam,b) for b in B])
bb=mul(C[3],(4,0))
cc=mul(add(C[2],neg(mul(bb,bb))),(4,0))
assert C==pm([cc,bb,one],[cc,bb,one])
# Exact target comparisons; all bounds are center densities.
from fractions import Fraction as F
p48=F(3,2)**24
cap48=F(10,4)**24 /5**24
assert cap48/p48==F(1,3)**24
bcs96=3**48*F(1,2)**24
cap96=F(14,4)**48 *2**24/7**48
assert cap96/bcs96==F(1,3)**48
print('PASS two density obstructions and literal short section.')

# Other CM twists: keep both elliptic Hodge lines even when their
# characters restrict to the same character of mu_2.
for p,q,k,want in [(5,25,2,{1:6,3:6}),(5,25,3,{1:2,3:4,5:2}),(7,49,2,{1:12,3:12})]:
    M=k*(q-1)
    hol={k*i-j for j in range(1,k) for i in range(1,(j*q-1)//k+1)}
    A={a for a in range(M) if a%k}
    assert all(2*sum(x in hol for x in O)==len(O) for O in orbits(A,p,M))
    old={(a,b) for a in A for b in (-1,1) if (a-b)%k==0}
    n=len(old)
    N=Counter()
    while old:
        a,b=min(old); u,v=a,b; s=[0]
        while (u,v) in old:
            old.remove((u,v))
            s.append(s[-1]+int(u in hol)+int(v==-1)-1)
            u,v=p*u%M,-v
        assert (u,v)==(a,b) and s[-1]==0
        N[sum(x-min(s) for x in s[:-1])]+=1
    assert n==2*(q-1) and N==want
    assert sum(a*z for a,z in N.items())==q-1
# Rational ratios from the three additional short morphisms.
assert F(10,4)**24*2**24/5**24/p48==F(2,3)**24
assert F(40,12)**24*3**12/5**24/p48==F(16,27)**12
assert F(14,4)**48*2**48/7**48/bcs96==F(2,9)**24
# Indices via elementary ramification counts of the explicit base quotients.
assert 2*2-2==2*(2*1-2)+2  # g(C)=2, involution has 2 fixed points
assert 2*4-2==4*(2*1-2)+3*2 # g(C)=4, three nontrivial involutions each fixes 2
assert 2*3-2==2*(2*1-2)+4 # g(C)=3, involution fixes 4
print('PASS the three other CM twists and their exact comparison arithmetic.')
