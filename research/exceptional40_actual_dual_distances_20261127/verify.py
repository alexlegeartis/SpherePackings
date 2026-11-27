#!/usr/bin/env python3
"""Exact four-cusp and finite quadratic checks, Python standard library."""
from fractions import Fraction as F
from itertools import product
from math import comb


def need(a):
    if not a:
        raise ArithmeticError('certificate check failed')


N, W = 11, 20  # retain Q^0,...,Q^10


def binom(n,k):
    return comb(n,k) if 0<=k<=n else 0


def mul(a,b):
    return [sum(a[i]*b[j-i] for i in range(j+1)) for j in range(N)]


def powers(a):
    out=[[1]+[0]*(N-1)]
    for _ in range(W):
        out.append(mul(out[-1],a))
    return out


x,y=[0]*N,[0]*N
for a,b in product(range(-8,9),repeat=2):
    n=3*(a*a+a*b+b*b)
    if n<N:
        x[n]+=1
    n+=3*(a+b)+1
    if n<N:
        y[n]+=1
need(x[:7]==[1,0,0,6,0,0,0]); need(y[:7]==[0,3,0,0,3,0,0])
X,Y=powers(x),powers(y)
P=[mul(X[W-j],Y[j]) for j in range(W+1)]
H=[[sum(binom(W-j,i-k)*2**(i-k)*binom(j,k)*(-1)**k
        for k in range(i+1)) for j in range(W+1)] for i in range(N)]


def solve(A,v):
    T=[[F(s) for s in a]+[F(b)] for a,b in zip(A,v)]
    n=len(v)
    for i in range(n):
        m=next(j for j in range(i,n) if T[j][i])
        T[i],T[m]=T[m],T[i]
        p=T[i][i]
        T[i]=[u/p for u in T[i]]
        for j in range(n):
            if j!=i and T[j][i]:
                s=T[j][i]
                T[j]=[a-s*b for a,b in zip(T[j],T[i])]
    return [r[-1] for r in T]


def jets(c,A5,n0,n1,n2):
    a=[F(0)]*(W+1)
    given=list(map(F,[c,0,0,0,0,A5]))
    for j in range(6):
        a[j]=(given[j]-sum(a[i]*P[i][j] for i in range(j)))/3**j
    for r,nr in enumerate([n0,n1,n2]):
        low=list(range(r,6,3))
        high=list(range(6+r,W+1,3))
        rhs=[729*nr*[1,0,0,F(-40,9),0][i]-sum(a[j]*H[i][j] for j in low)
             for i in range(5)]
        z=solve([[H[i][j] for j in high] for i in range(5)],rhs)
        for j,b in zip(high,z):
            a[j]=b
    f=[sum(a[j]*P[j][i] for j in range(i+1)) for i in range(N)]
    b=[sum(a[j]*H[i][j] for j in range(W+1)) for i in range(N)]
    g=[sum(b[j]*P[j][i] for j in range(i+1)) for i in range(N)]
    need(f[:6]==given)
    need(g[:5]==[729*(n0+n1+n2),0,0,0,0])
    return [f[6],f[7],f[8],f[9],g[5]]


columns=[jets(*(int(i==j) for j in range(5))) for i in range(5)]
need(columns==[
    [-495,0,0,7410,22320522],
    [0,0,-40,0,59049],
    [495,0,0,1084590,-275562],
    [0,9180,0,0,551124],
    [0,0,116235,0,-1102248],
])
need(jets(1,0,1,0,0)[3:]==[1092000,729*30240])

# Isotropic slice of E8/3 E8 in the doubled-coordinate presentation.
u=[(1,1,1,0,0,0,0,0),(1,2,0,1,0,0,0,0),
   (0,0,0,0,1,1,1,0),(0,0,0,0,1,2,0,1)]
def dot(x,y): return sum(a*b for a,b in zip(x,y))%3
need(all(dot(a,b)==0 for a,b in product(u,repeat=2)))
sub={tuple(sum(c*v[i] for c,v in zip(w,u))%3 for i in range(8))
     for w in product(range(3),repeat=4)}
need(len(sub)==81 and 3**8//len(sub)==81)

# Gram I_m+J_m for the constant q_K=1 finite distance case.
# If m>=9, m=9 already has rank 9 in F3 and is impossible.
need((1+9)%3!=0)
for mu,m in [(F(16,3),9),(F(14,3),81),(F(4),2187),(F(10,3),6561)]:
    need(m*(mu/6)**20<1)
need(F(3,2)**20/81==F(3**16,2**20))

# Degree-one and degree-two coset harmonic factors have only j=0
# at T-residue 2 after dividing by the universal order-5 cusp factor.
need(all([j for j in range(d+1) if j%3==0]==[0] for d in [1,2]))
print('PASS: 40D coset theta jets, density barriers, isotropic cap, deep-phase harmonic degree check.')
