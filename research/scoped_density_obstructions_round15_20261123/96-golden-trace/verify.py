"""Exact scalar identities for the golden trace norm and endpoints.

Run python -B verify.py. No third-party packages.
"""
from fractions import Fraction as F
from math import prod

N = 7

def mul(f, g):
    return [sum(f[j]*g[i-j] for j in range(i+1)) for i in range(N+1)]

def power(f, k):
    if not k:
        return [1]+[0]*N
    return mul(f,power(f,k-1))

def sigma(n,k):
    return sum(d**k for d in range(1,n+1) if n%d==0)

def subtract(f,g,c):
    return [a-c*b for a,b in zip(f,g)]

E4 = [1]+[240*sigma(i,3) for i in range(1,N+1)]
E6 = [1]+[-504*sigma(i,5) for i in range(1,N+1)]
Delta = [F(c,1728) for c in subtract(power(E4,3),power(E6,2),1)]
basis = [mul(power(Delta,j),power(E4,12-3*j)) for j in range(5)]
ext = basis[0]
for j in range(1,5):
    ext = subtract(ext,basis[j],ext[j])
assert ext[:6] == [1,0,0,0,0,565866362880]
assert basis[-1][4:6] == [1,-96]

def harmonic_polynomial(k):
    # Coefficient list for q^(2i) s^(k-2i), leading coefficient 1.
    p=[F(1)]
    for i in range(1,k//2+1):
        a,b=2*i,k-2*i
        p.append(-p[-1]*F(20*(b+2)*(b+1),2*a*(2*a+94+4*b)))
    for i in range(1,len(p)):
        a,b=2*i,k-2*i
        assert p[i] * (2*a*(2*a+94+4*b)) + p[i-1]*20*(b+2)*(b+1) == 0
    return p

def harmonic(k,q,s):
    return sum(p*q**(2*i)*s**(k-2*i) for i,p in enumerate(harmonic_polynomial(k)))

assert harmonic_polynomial(2)==[1,-F(5,49)]
assert harmonic_polynomial(4)==[1,-F(30,53),F(25,901)]
assert [mul(basis[-1],power(E4,j))[5] for j in (1,2)] == [144,384]
c = [462780864000, 99773856000, 3311642880]
a = [-40,-60,4]
scls = [2,6,10]
assert sum(c)==ext[5] and sum(a)==-96
for k,modfactor in ((2,144),(4,384)):
    assert sum(harmonic(k,10,s)*b for s,b in zip(scls,c))==0
    assert sum(harmonic(k,10,s)*b for s,b in zip(scls,a))==modfactor*harmonic(k,8,0)
    assert len({harmonic(k,10,s) for s in scls})==3
# Upper bound on A, confirmed by nonnegativity of middle type.
assert F(c[1],-a[1])==1662897600
# Odd H_1,H_3,H_5 are triangular in s,s^3,s^5, Vandermonde nonzero.
assert prod(scls)*prod(y*y-x*x for i,x in enumerate(scls) for y in scls[i+1:])!=0

# Pairing grids follow s congruent to q modulo 4 and min bounds for U, U^-1.
for m,exp in ((4,[-4,0,4]),(6,[-6,-2,2,6]),(8,[-8,-4,0,4,8]),(10,[-10,-6,-2,2,6,10])):
    assert [s for s in range(-m,m+1) if (s-m)%4==0]==exp
    assert all(((3*m+s)//2)>=m and ((3*m-s)//2)>=m for s in exp)
# If all norm-8 pairings vanish, s= +/-14 at q=10 would give forbidden norm-8 pairings.
assert ((30-14)//2,(50-3*14)//2)==(8,4)
assert ((30-14)//2,(-50+3*14)//2)==(8,-4)
# Endpoint q + s/5 = 8 gives only these possibilities for q>=8.
assert [(q,40-5*q) for q in range(8,42,2) if abs(40-5*q)<=3*q-16]==[(8,0),(10,-10),(12,-20)]
assert ((36-20)//2,(60-60)//2)==(8,0)
# For extremal parents in dimensions 24/48/72, variance or kurtosis forces |s|=m.
assert F(5*4**2,13)>0
assert F(5*6**2,25)>2**2
assert F(3*25*8**4,37*39)>16*F(5*8**2,37)
# Both signs for |s|=m follow the units, even without a Galois lattice automorphism.
for m in (4,6,8,10):
    assert ((3*m-m)//2,(5*m-3*m)//2)==(m,m)
# Squared Hermite maxima of min-8 classes in order h=0,4,8.
assert [max(F(64),F((40-h)**2,20)) for h in (0,4,8)] == [80,F(324,5),64]
# In the min-8 canceled branch, endpoint kissing count is 2A + c[-1]/2 + a[-1] A/2.
assert c[-1]//2 == 1655821440 and 2+a[-1]//2 ==4
print("PASS: golden trace endomorphism and shell arithmetic")

# Rank 88 transfer: rebuild the ordinary coefficients and harmonic ratios.
bb=[mul(power(Delta,j),power(E4,11-3*j)) for j in range(4)]
tt=bb[0]
for j in range(1,4):
    tt=subtract(tt,bb[j],tt[j])
assert tt[:6]==[1,0,0,0,168498000,2480127344640]
assert [mul(power(Delta,4),power(E4,j))[5] for j in (0,1)]==[-96,144]
def H88(k,q,s):
    p,res=F(1),F(s**k)
    for i in range(1,k//2+1):
        a,b=2*i,k-2*i
        p *= -F(20*(b+2)*(b+1),2*a*(2*a+86+4*b))
        res += p*q**a*s**b
    return res
count88=[1970971200000,486367200000,22788944640]
assert sum(count88)==tt[5]
for k,ratio in ((2,-96),(4,144)):
    assert sum(a*H88(k,10,s) for a,s in zip(count88,scls))==tt[4]*ratio*H88(k,8,0)
assert 2*tt[4]+count88[-1]//2==11731468320
assert F(180,5) < (F(5,4)**22*F('4.95809518917e-20')/F('1.85638056287e-19')) < 37
assert F(27) < (F(5,4)**22*F('4.95809518917e-20')/F('2.42250297455694e-19')) < 28
print("PASS: rank-88 golden transfer shell arithmetic")

# Rank 96: all absolute pairing types in the next shell.
S12=(0,4,8,12,16,20)
C12=(1218201435360000,1453085210304000,286338229632000,14383631808000,99773856000,0)
A12=(3212,-1540,2480,370,-60,2)
assert sum(C12)==ext[6]==2972108280960000
assert sum(A12)==basis[-1][6]==4464
for k in (2,4,6):
    series=mul(power(Delta,4),power(E4,k//2))
    for A,Vs,Ps in ((0,C12,c),(1,A12,a)):
        M4=A*harmonic(k,8,0)
        M5=sum(w*harmonic(k,10,s) for w,s in zip(Ps,scls))
        M6=M4*series[6]
        if k==6:
            M6+=(M5-M4*series[5])*power(Delta,5)[6]
        else:
            assert M5==M4*series[5]
        assert sum(b*harmonic(k,12,z) for b,z in zip(Vs,S12))==M6
assert prod(b*b-a*a for i,a in enumerate(S12[:4]) for b in S12[i+1:4])
assert (C12[-2],A12[-2])==(c[1],a[1])
assert ((36-16)//2,(60-3*16)//2)==(10,6)
print("PASS: rank-96 norm-12 pairing split")
