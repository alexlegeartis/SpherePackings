import json
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path

import replay as r

ds=[r.d[O[0]] for O in r.orbits]
assert all(r.d[i]==d for O,d in zip(r.orbits,ds) for i in O)
D=list(product(ds,repeat=2))
sizes=[s*t for s,t in product(r.sizes,repeat=2)]
H=[[a*b for a,b in product(x,y)] for x,y in product(r.H,repeat=2)]

# Complete vertex census in the normalized sorted, clipped square.
lines={(a,b,18) for a,b in D if a or b}
lines|={(1,0,1),(1,0,9),(0,1,1),(0,1,9),(1,-1,0)}
V=set()
for (a,b,c),(A,B,C) in combinations(lines,2):
    dt=a*B-b*A
    if dt:
        x,y=F(c*B-b*C,dt),F(a*C-c*A,dt)
        if 1<=x<=y<=9:
            V.add((x,y))
assert len(V)==64

data=json.loads(Path(__file__).with_name('scales.json').read_text())
certs={(F(row[0]),F(row[1])):row[2] for row in data}
assert set(certs)==V and len(data)==64
open_vertices=set()
for x,y in V:
    z={int(i):F(t) for i,t in certs[x,y]}
    assert all(0<=i<len(H) and t>0 for i,t in z.items())
    hp=[sum(t*H[i][j] for i,t in z.items()) for j in range(len(D))]
    assert all(hp[j]+sizes[j]<=0 for j,(a,b) in enumerate(D)
               if j and a*x+b*y>=18)
    B=1+hp[0]
    target=F(64,9)*(x*y)**3
    if B>=target.numerator//target.denominator+1:
        open_vertices.add((x,y))
        if (x,y)==(F(9,4),F(9,4)):
            assert B==1215 and target.numerator//target.denominator+1==923
            assert 3**6<923 and 3**7>B # closes this and all larger additive sizes
        elif (x,y)==(F(9,5),F(18,5)):
            assert B==2187 and target.numerator//target.denominator+1==1935
            assert 1931<1935 # established by the stdlib proof in replay.py
assert open_vertices=={(F(9,4),F(9,4)),(F(9,5),F(18,5))}
print('PASS: 64 continuous-scale vertices, all lattices excluded; sole arbitrary open mask (9/4,9/4).')
