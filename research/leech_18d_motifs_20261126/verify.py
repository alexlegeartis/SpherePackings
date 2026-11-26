#!/usr/bin/env python3
"""Exact baseline replay for two primitive rank-six Leech embeddings."""
from itertools import product, combinations, permutations
from collections import Counter
from math import prod
from fractions import Fraction as F
from pathlib import Path
import json

bases=json.loads(Path(__file__).with_name('bases.json').read_text())
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def mv(A,v):return tuple(dot(row,v) for row in A)
def det(A):
    n=len(A)
    return sum((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
               *prod(A[i][p[i]] for i in range(n)) for p in permutations(range(n)))
def adj(A):
    n=len(A)
    return [[(-1)**(i+j)*det([[A[k][l] for l in range(n) if l!=i]
                              for k in range(n) if k!=j])
             for j in range(n)] for i in range(n)]

def golay():
    for t in range(4096):
        c=0
        for i in range(12):
            if t>>i&1:c^=0xae3<<i
        yield c|((c.bit_count()%2)<<23)

def shell():
    C=list(golay())
    assert len(set(C))==4096
    assert Counter(x.bit_count() for x in C)=={0:1,8:759,12:2576,16:759,24:1}
    for i,j in combinations(range(24),2):
        for a,b in product((-4,4),repeat=2):
            v=[0]*24;v[i],v[j]=a,b
            yield v
    for c in C:
        if c.bit_count()==8:
            I=[i for i in range(24) if c>>i&1]
            for s in range(256):
                if s.bit_count()%2:continue
                v=[0]*24
                for j,i in enumerate(I):v[i]=-2 if s>>j&1 else 2
                yield v
    for c in C:
        base=[-1 if c>>i&1 else 1 for i in range(24)]
        for i in range(24):
            v=base[:];v[i]*=-3
            yield v

qhist=[{0:1,162:72,324:90,486:80},{0:1,128:54,192:36,320:72,384:27,512:2}]
best=[]
for b,B in enumerate(bases):
    G=[[dot(a,c)//8 for c in B] for a in B]
    assert all(dot(a,c)%8==0 for a in B for c in B)
    assert all(G[i][i]==4 for i in range(6))
    D=det(G);J=adj(G)
    assert D==[243,192][b]
    R=[[4,-2,0,0,0,0],[-2,4,-2,0,0,0],[0,-2,4,-2,0,-2],
       [0,0,-2,4,-2,0],[0,0,0,-2,4,0],[0,0,-2,0,0,4]]
    if b==0:
        A=[[0,0,1,0,-1,-1],[0,0,0,1,-1,-1],[0,0,1,-1,0,-1],
           [1,0,0,0,0,-1],[0,1,0,0,0,-1],[1,-1,1,0,0,-1]]
        C=adj([[t//2 for t in row] for row in R])
        assert abs(det(A))==1 and [[dot(x,mv(C,y)) for y in A] for x in A]==G
    else:assert G==R
    assert all(sum(J[i][k]*G[k][j] for k in range(6))==(D if i==j else 0)
               for i in range(6) for j in range(6))
    key=lambda j:tuple(a%D for a in j)
    zero=(0,)*6
    q={}
    for z in product(range(-3,4), repeat=6):
        j=mv(J,z);k=key(j);u=dot(j,z)
        if u<q.get(k,100000):q[k]=u
    # Cauchy |z_i|^2 <= 4 q: if shortest found <= 8/3, no omitted
    # representative (one |z_i| >= 4) can improve any class.
    assert len(q)==D and max(q.values())<=F(8*D,3)
    assert Counter(q.values())==qhist[b]
    atmin=set(); members=set(); size=0
    for y in shell():
        if y in B:members.add(tuple(y))
        tmp=mv(B,y);assert all(x%8==0 for x in tmp)
        z=tuple(a//8 for a in tmp);j=mv(J,z);k=key(j)
        if dot(z,j)==q[k]:atmin.add(k)
        size+=1
    assert size==196560 and len(members)==6 and len(atmin)==D
    def delta(a,c):return q[tuple((x-y)%D for x,y in zip(a,c))]
    def maxclique(u):
        P=[x for x,t in q.items() if 0<t<=u]
        adj={x:{y for y in P if x!=y and delta(x,y)<=u} for x in P}
        def go(M,P):
            n=M
            while P and M+len(P)>n:
                x=P.pop()
                n=max(n,go(M+1,P&adj[x]))
            return n
        return go(1,set(P))
    if b==0:
        assert all(all(u%81==0 for u in p) for p in q)
        packet=((0,0,0,0,0,0),(0,1,1,1,0,0),(-1,0,0,1,0,0),
                (0,1,0,1,1,0),(0,0,1,1,1,0),(-1,0,0,0,0,-1))
        keys=[key(mv(J,x)) for x in packet]
        assert len(set(keys))==6 and all(delta(a,c)==162 for a,c in combinations(keys,2))
        assert maxclique(162)==6
        # The 80 forbidden isotropic differences, and all restrictions of
        # ternary coordinate characters, give an exact Hoffman bound 27.
        isotropic=[tuple(u//81 for u in x) for x in q if q[x]==486]
        E=set()
        for y in product(range(3),repeat=6):
            c=Counter(dot(x,y)%3 for x in isotropic);assert c[1]==c[2]
            E.add(c[0]-c[1])
        assert E=={-10,-1,8,80}
        bound2=[F(1,D),F(36,D)*F(5,6)**18,F(27**2,D)*F(2,3)**18,
                D*F(1,2)**18]
        assert max(bound2)==bound2[1]
    else:
        assert maxclique(128)==3 and maxclique(192)==9
        sector=lambda p:tuple(2*u%D for u in p)
        assert len({sector(p) for p in q})==3
        assert all(q[p] in (0,192,384) for p in q if sector(p)==zero)
        bound2=[F(1,D),F(9,D)*F(5,6)**18,F(81,D)*F(3,4)**18,
                F(24**2,D)*F(7,12)**18,D*F(1,2)**18]
        assert max(bound2)==bound2[0]
    assert max(bound2)<F(3,4)**18
    best.append(float(max(bound2))**.5)
    print('det',D,'classes',len(q),'exact norm-four at minima',len(atmin),
          'spectrum',sorted(Counter(q.values()).items()),'center cap',best[-1])
print('PASS both exact minima, primitivity, clique/radius ceilings, and six-coset witness')
