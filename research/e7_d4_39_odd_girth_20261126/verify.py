#!/usr/bin/env python3
"""Independent exact 39D radial-alphabet/odd-girth replay; standard library."""
from itertools import combinations
from math import isqrt
from pathlib import Path
from collections import Counter
import json

cert=json.loads(Path(__file__).with_name('e7_input.json').read_text())
G=cert['E7_gram']
M=cert['phase_isometry']
alpha=cert['alpha']
assert cert['field_polynomial']==67 and len(set(alpha))==5
Di=[120,90,80,75,48,45,40]
Li=[[60,0,0,0,0,0,0],[-30,60,0,0,0,0,0],[-30,-20,60,0,0,0,0],
    [0,0,-45,60,0,0,0],[-30,-20,-30,-24,60,0,0],
    [0,0,0,0,-75,60,0],[0,0,0,0,0,-80,60]]
A=216000
assert all(sum(Li[i][k]*Li[j][k]*Di[k] for k in range(7))==A*G[i][j]
           for i in range(7) for j in range(7))

T=[255]*4**7
LIM=40*A
z=[0]*7
def go(j,rem,p):
    if j<0:
        assert (LIM-rem)%A==0
        T[p]=min(T[p],(LIM-rem)//A)
        return
    w=sum(Li[i][j]*z[i] for i in range(j+1,7))
    r=isqrt(rem//Di[j])
    for a in range(-((r+w)//60),(r-w)//60+1):
        z[j]=a
        go(j-1,rem-Di[j]*(60*a+w)**2,p|((a&3)<<(2*j)))
go(6,LIM,0)
assert max(T)==24

def pack(f):return sum(1<<(2*j) for j in range(7) if f>>j&1)
assert sorted(M)==list(range(128)) and M[64]==82 and M[0]==0
rcoords=[(82>>j)&1 for j in range(7)]
dot=[sum(x*y for x,y in zip(row,rcoords)) for row in G]
assert all(v%2==0 for v in dot) and sum(x*y for x,y in zip(dot,rcoords))==6
assert 5*128**2*2**40 < 3**32*9**2  # 128/9 < pinned 39D incumbent
assert 5*8192**2*2**40 > 3**32*441**2 # 8192/441 would beat it
C=[[T[pack(f)^(2*pack(M[g]))] for g in range(128)] for f in range(128)]
projected=[[min(C[f][u],C[f][u+64]) for u in range(64)] for f in range(128)]
types=[min(c) for c in projected]
assert Counter(types)=={0:1,2:63,4:63,6:1} and types[82]==6
rad=[u for u in range(64) if projected[82][u]==6]
assert len(rad)==28
assert all(C[82][u]==C[82][u^64] and C[82][u]==(6 if u in rad else 14) for u in range(128) if u<64)
roots=[f for f in range(128) if types[f]==2]
evens=[f for f in range(128) if types[f]==4]
assert all(sum(c==2 for c in projected[f])==2 and sum(c<=6 for c in projected[f])==34 for f in roots)
assert all(sum(c==4 for c in projected[f])==12 and sum(c<=8 for c in projected[f])==44 for f in evens)

def mul(x,y):
    v=0
    while y:
        if y&1:v^=x
        x<<=1
        if x&64:x^=67
        y>>=1
    return v
def q(x):
    y=1
    for _ in range(9): y=mul(y,x)
    yy=mul(y,y)
    return y^yy^mul(yy,yy)
assert set(rad)=={x for x in range(64) if q(x^60)==0}
assert all(M[x^y]==M[x]^M[y] for x in range(128) for y in range(128))
assert all(any(len({(x&d).bit_count()%2 for x in range(64)
                         if projected[f][x]==6})==1 for d in range(1,64)) for f in roots)
A0=[[(u | (mul(a,u)<<6)) for u in range(64)] for a in alpha]

for i,j in combinations(range(5),2):
    a,b,c=set(range(5))-{i,j}
    covered={A0[a][x]^A0[b][y]^A0[c][z] for x in rad for y in rad for z in rad}
    holes=set(range(4096))-covered
    assert len(holes)==36
    lines=[]
    for k in (i,j):
        # A nonempty displacement of any hole set into an insufficient line
        # would take some source point into a hole; anchor at the least one.
        mx=max(Counter((s>>6)^mul(alpha[k],s&63) for s in holes).values())
        lines.append(mx)
        for f in roots:
            cheap=[A0[k][u] for u in range(64) if projected[f][u]<=6]
            assert all(any(s^cheap[0]^p in covered for p in cheap) for s in holes)
    assert lines in ([1,1],[36,1],[1,36]),lines
    assert len({v^A0[i][x] for v in covered for x in rad})==4096
    assert len({v^A0[j][x] for v in covered for x in rad})==4096
    print('radical triple slopes',sorted(alpha[k] for k in (a,b,c)),
          'holes',len(holes),'max intersections with omitted slopes',lines)

# Explicit 15-cycle in the full half-grid quotient. All high edges are at a
# shift with four r/2 blocks, and one local even nonzero half-grid block.
r=82
u,v,w=next((a,b,a^b^r) for a,b in combinations(roots,2) if a^b^r in roots)
assert u^v^w==r and len({u,v,w})==3
phase_edges=[]
for block in range(5):
    for t in (u,v,w):
        row=[r]*5;row[block]=r^t
        assert types[row[block]]==4
        phase_edges.append(row)
state=(0,)*5
seen={state}
for k,e in enumerate(phase_edges,1):
    state=tuple(a^b for a,b in zip(state,e))
    if k<15:
        assert state not in seen
        seen.add(state)
    else:assert state==(0,)*5
# Coefficient numerators modulo 4 can be set arbitrarily on the last edge
# subject to the prescribed fractional mask. Choose its opposite total.
sums=[[sum(e[b]>>j & 1 for e in phase_edges[:-1]) for j in range(7)] for b in range(5)]
last=[[(-s)%4 for s in row] for row in sums]
assert [sum((k%2)<<j for j,k in enumerate(row)) for row in last]==phase_edges[-1]
print('witness 15-cycle deviations',u,v,w,'all 15 projected vertices distinct')
print('PASS local minima, 10 three-bank covers, every four-bank cover, and odd-girth witness')

# Exact 11-cycle one level beyond the half-grid. The basis below has
# twice the coordinates of A7 plus the half-integral sum-zero coset.
def ai(i):
    return [2*((j==i)-(j==i+1)) for j in range(8)]
b2=[ai(2),[-1,-1,-1,1,1,1,1,-1],ai(1),ai(0),ai(3),ai(4),ai(5)]
assert all(sum(x*y for x,y in zip(b2[i],b2[j]))==4*G[i][j] for i in range(7) for j in range(7))
rvec=[(82>>j)&1 for j in range(7)]
w4=[4,1,4,4,3,2,1]
assert [sum(w4[i]*b2[i][j] for i in range(7)) for j in range(8)]==[7]+[-1]*7
assert all((2*a-2*b)%4==0 for a,b in zip(w4,rvec)) # 2w = R mod E7
# Numerators in eighth-units. Closest independent rounding to the integer
# or half-integer coordinate classes; enforce sum zero.
target=[7]+[-1]*7
assert sum(b*b for b in target)==56
for shift,goal in [(0,56),(4,120)]:
    nearest=[min(range(shift-16,shift+24,8),key=lambda k:abs(k-z)) for z in target]
    err=[z-y for z,y in zip(target,nearest)]
    assert len(set(err))==1 and sum(nearest)%8==0
    d=sum(z*z for z in err)
    n=-(sum(nearest)//8)
    # The quadratic penalties for further steps at same coordinate
    # increase, so distinct coordinates give the optimum for |n|<=8.
    step=8 if n>0 else -8
    assert d+abs(n)*((err[0]-step)**2-err[0]**2)==goal
r4=[2*x for x in rvec]
edges=[]
for b in range(5):
    row=[r4]*5;row[b]=w4
    edges.extend([row,row])
prior=[[sum(e[i][j] for e in edges) for j in range(7)] for i in range(5)]
last=[[(-v)%8 for v in z] for z in prior]
assert all([x%4 for x in z]==[x%4 for x in r4] for z in last)
assert all((a+b)%8==0 for x,y in zip(prior,last) for a,b in zip(x,y))
edges.append(last)
state=[[0]*7 for _ in range(5)]
seen={tuple(map(tuple,state))}
for k,e in enumerate(edges,1):
    state=[[(u+v)%4 for u,v in zip(x,y)] for x,y in zip(state,e)]
    flat=tuple(map(tuple,state))
    if k<11:
        assert flat not in seen
        seen.add(flat)
    else:assert not any(map(any,state))
print('witness 11-cycle at squared threshold 55/8: 4w =',w4)
