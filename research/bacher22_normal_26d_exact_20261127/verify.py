#!/usr/bin/env python3
"""Exact arithmetic replay for the Bacher-22 projection bounds."""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb


def need(a):
    if not a:
        raise ArithmeticError('certificate check failed')


W = [0,3,4,3,2,3,2,7,8,7,4,7,2,7,2,3,
     4,3,8,3,2,7,2,3,4,7,8,7,2,3,2,7,
     4,3,8,3,2,7,2,3,4,7,8,7,2,3,2,7,
     4,7,8,7,6,7,6,3,4,3,8,3,6,3,6,7]


def add(a,b):
    return ((a^b)^(((a&b)&5)<<1))&63


def neg(a):
    return add(a,add(a,a))


def parity(a):
    return (a&1)|((a&4)>>1)|((a&48)>>2)


short = [a for a,w in enumerate(W) if w==2]
need({parity(a) for a in short}=={2,6,10})
for a,b in combinations(short,2):
    if parity(a)!=parity(b):
        need(W[add(a,b)]>=4 and W[add(a,neg(b))]>=4)

# The short graph DOES have four-cliques.  None also satisfies the
# three opposite-edge sum/difference bounds for an A3 packet.
close = [a for a,w in enumerate(W) if 0<w<4]
tetra = [x for x in combinations(close,3)
         if all(0<W[add(a,neg(b))]<4 for a,b in combinations(x,2))]
need(len(tetra)==144 and (1,5,6) in tetra)
need(not any(all(0<W[add(a,neg(b))]<4 for a,b in combinations(x,2))
             for x in combinations(close,4)))
for a,b,c in tetra:
    opposites = [add(add(a,b),neg(c)),add(add(a,c),neg(b)),add(add(b,c),neg(a))]
    need(any(x==0 or W[x]>=7 for x in opposites))

# Each polynomial below is in z=4 lambda, with 0<=z<=1.
N = 27


def p(a=0,b=0):
    return [F(a), F(b)] + [F(0)]*(N-2)


def plus(a,b):
    return [x+y for x,y in zip(a,b)]


def scale(a,b):
    return [x*b for x in a]


def mul(a,b):
    return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(N)]


def mul_all(xs):
    out=p(1)
    for a in xs:
        out=mul(out,a)
    return out


def det_poly(g):
    r=len(g)
    out=p()
    for w in permutations(range(r)):
        sign=(-1)**sum(w[i]>w[j] for i in range(r) for j in range(i+1,r))
        out=plus(out,scale(mul_all(g[i][w[i]] for i in range(r)),sign))
    return out


lam=p(0,F(1,4))
half=p(F(1,2))
upper=p(F(1,4))
zero=p()
classes=[
    [[]],
    [[zero],[upper]],
    [[zero,lam],[zero,upper],[upper,upper]],
    [[zero,lam,lam],[zero,lam,upper],[zero,upper,upper],
     [upper,upper,upper],[scale(lam,F(1,2))]*3],
]
inverse_H=[F(1,4),F(1,2),F(3,4),F(1)]
against=scale(mul_all([p(1,F(-1,16))]*26),F(3,64))
zeros=[]
for r,vertices in enumerate(classes):
    for i,x in enumerate(vertices):
        a=[plus(plus(half,scale(lam,-1)),t) for t in x]
        f=p(F(3,4),F(-1,4))
        for v in a:
            f=plus(f,scale(v,F(-1,4)))
        g=[[a[j] if j==k else scale(plus(plus(x[j],x[k]),scale(lam,-1)),F(-1,2))
            for k in range(r)] for j in range(r)]
        p1=scale(mul(det_poly(g),mul_all([f]*(4-r))),inverse_H[r])
        d=plus(p1,scale(against,-1))
        b=[sum(d[j]*F(comb(k,j),comb(26,j)) for j in range(k+1)) for k in range(27)]
        need(min(b)>=0 and min(b[1:])>0)
        if b[0]==0:
            zeros.append((r,i))
need(zeros==[(2,0),(3,0),(3,3),(3,4)])

# Low-radius 26D branch, and all-periodic 28D reduction using the
# published 6D center-density upper bound rounded upward to 0.081.
need(3*F(15,16)**22<1)
need(8*F(81,1000)<F(2,3))
need((8*F(1,8))**2/F(2)==F(1,2))  # 5D D5 center density 1/(8 sqrt(2))
need(4*F(7,16)**5>F(3,64) and 6*F(7,16)**5>F(3,64))

# Exact failure of one proposed 27D oblique-tail placement: coordinates
# (0,0,0,-2,1) have zero label and length^2 7/2 rather than at least 4.
g=[4,20,36,1,2]
need(add(neg(add(g[3],g[3])),g[4])==0)
need(F(3,8)*4+1-4*F(1,8)+F(1,2)*3==F(7,2))
print('PASS: finite projected-label identities, all 11 Bernstein certificates, low-radius/periodic comparisons, oblique-tail counterexample.')
