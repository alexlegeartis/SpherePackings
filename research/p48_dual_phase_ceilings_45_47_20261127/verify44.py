#!/usr/bin/env python3
"""Exact membership and clique-bound replay for the literal det-605 normal."""
from itertools import combinations, product
from pathlib import Path
from fractions import Fraction as F
import json

W=json.loads((Path(__file__).parent/"t44_witnesses.json").read_text())
B=W["normals"]
assert W["scale"]==12 and len(B)==4 and all(len(x)==48 for x in B)
G=[[6,-2,-1,-1],[-2,6,-2,-1],[-1,-2,6,-2],[-1,-1,-2,6]]
M=[[14,8,7,6],[8,16,9,7],[7,9,16,8],[6,7,8,14]]
assert all(sum(x*y for x,y in zip(a,b))==12*G[i][j]
           for i,a in enumerate(B) for j,b in enumerate(B))
assert all(sum(G[i][k]*M[k][j] for k in range(4))==55*(i==j)
           for i,j in product(range(4),repeat=2))
def det(A):
    if len(A)==1:
        return A[0][0]
    return sum((-1)**j*x*det([r[:j]+r[j+1:] for r in A[1:]]) for j,x in enumerate(A[0]))
assert det(G)==605

QR={j*j%23 for j in range(1,23)}
S=[[0 if i==j else 1 if i==0 else -1 if j==0 else 1 if (i-j)%23 in QR else -1
    for j in range(24)] for i in range(24)]
def code(x):
    return all((x[24+j]-sum(S[i][j]*x[i] for i in range(24)))%3==0
               for j in range(24))
def member(y):
    parity=y[0]%2
    if any(x%2!=parity for x in y):
        return False
    x=[(a-3*parity)//2 for a in y]
    return sum(x)%2==parity and code(x)
assert code([-2]*48) and all(member(y) for y in B)
def key(a):
    return tuple(sum(x*y for x,y in zip(r,a))%55 for r in M)
def norm(a):
    return sum(x*sum(y*z for y,z in zip(r,a)) for x,r in zip(a,M))

C={}
for a in product(range(-6,7),repeat=4):
    q,k=norm(a),key(a)
    C[k]=min(q,C.get(k,10**6))
assert len(C)==605 and max(C.values())==176
# All excluded labels have a coordinate >=7 in magnitude, forcing normal norm >=49/6.
assert F(49,6)>F(176,55)
D={}
for n,count in [(6,547),(8,58)]:
    before=len(D)
    for s in W["norm"+str(n)]:
        y=[ord(v)-ord("g") for v in s]
        assert len(y)==48 and member(y) and sum(x*x for x in y)==12*n
        for p in (y,[-x for x in y]):
            t=[sum(x*u for x,u in zip(b,p)) for b in B]
            assert all(a%12==0 for a in t)
            a=tuple(x//12 for x in t)
            k=key(a)
            assert C[k]==norm(a)
            assert k not in D or not any(k)
            D[k]=55*n-C[k]
    assert len(D)-before==count
assert D[0,0,0,0]==330
assert D.keys()==C.keys()

def diff(a,b):
    return tuple((x-y)%55 for x,y in zip(a,b))
def at_threshold(q,cap=None):
    goal=(5*314**22)//q**22 if cap is None else cap
    # An improving clique would have goal+1 points; translate one to zero.
    v=[k for k,n in D.items() if any(k) and n>=q]
    if len(v)<goal:
        return
    adj=[]
    for i,a in enumerate(v):
        adj.append(sum(1<<j for j,b in enumerate(v) if i!=j and D[diff(a,b)]>=q))
    def exists(todo,k):
        if k==0:
            return True
        if todo.bit_count()<k:
            return False
        remain=todo
        order=[]
        c=0
        while remain:
            available=remain
            c+=1
            while available:
                bit=available & -available
                i=bit.bit_length()-1
                order.append((i,c))
                remain^=bit
                available&=~(bit|adj[i])
        for i,c in reversed(order):
            if c<k:
                return False
            if exists(todo&adj[i],k-1):
                return True
            todo^=1<<i
        return False
    assert not exists((1<<len(v))-1,goal), (q,goal)

for q in set(D.values()):
    at_threshold(q)
at_threshold(316,2)
S=[(0,0,0,0)]+[tuple(int(i==j) for j in range(4)) for i in range(4)]
assert len({key(a) for a in S})==5
assert max(norm(tuple(x-y for x,y in zip(a,b))) for a,b in combinations(S,2))==16
assert F(25,605)*F(314,55*4)**44==F(157**44,2**44*5**43*11**46)
print("PASS: 605 explicit phase lifts and sharp scoped 44D density")

