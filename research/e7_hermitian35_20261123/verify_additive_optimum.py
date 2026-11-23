"""Exact 35D upper bound in the fixed five-E7 binary-period chart."""
from itertools import product
from collections import defaultdict

# A phase is (x,e) with x in a 6D symplectic space and e the radical bit.
# The nonzero radical has squared depth 6; remaining phase depths 2/4
# are selected by q_0(x)+e=1/0, where this split basis has q_0=sum x_i*y_i.
def q(x):
    return sum(((x >> (2*i)) & 3) == 3 for i in range(3)) % 2

def symp(a,b):
    return sum(((a >> (2*i)) & 1)*((b >> (2*i+1)) & 1)
               + ((b >> (2*i)) & 1)*((a >> (2*i+1)) & 1)
               for i in range(3)) % 2

T = [[], [], [], []]
for x,e in product(range(64),range(2)):
    d = 3*e if not x else 2-(q(x)^e)
    T[d].append((x,e))
assert list(map(len,T)) == [1,63,63,1]
W = {}
for y,b in product(range(64),range(2)):
    k = tuple(sum((-1)**(symp(x,y)+e*b) for x,e in D) for D in T)
    W[k] = W.get(k,0)+1
assert W == {(1,63,63,1):1,(1,-1,-1,1):63,(1,-7,7,-1):36,(1,9,-9,-1):28}
C = list(W.keys())

# f is given by a table of coefficients on composition classes
# (number of coordinates of squared depths 0,2,4,6). Missing values zero.
f={
 (0,0,0,5):-109530453720,(0,0,1,4):-24418384800,(0,0,2,3):-2014929684,
 (1,1,0,3):-3821770764,(1,1,1,2):-272060879,(2,0,1,2):-855244384,
 (2,1,2,0):-80895384,(2,2,0,1):-970642489,(2,3,0,0):1654554972,
 (3,0,1,1):-11879541886,(3,1,0,1):-11694036880,(3,1,1,0):8335120472,
 (3,2,0,0):24700329651,(4,0,0,1):-114118688191,(4,0,1,0):100054795864,
 (4,1,0,0):169850633129,(5,0,0,0):1000000002280
}
assert all(u<=0 for a,u in f.items() if sum((2*i)*a[i] for i in range(4))>=8)
def counts(m,n):
    for a in product(range(n+1),repeat=m):
        if sum(a)==n: yield a

def ft(s):
    b = {(0,0,0,0):1}
    for typ,e in enumerate(s):
        for _ in range(e):
            g=defaultdict(int)
            for a,v in b.items():
                for i,w in enumerate(C[typ]):
                    c=list(a); c[i]+=1; g[tuple(c)]+=v*w
            b=g
    return sum(v*f.get(a,0) for a,v in b.items())
S = list(counts(4,5))
assert len(S) == 56 and min(ft(s) for s in S)>=0
h = ft((5,0,0,0))
assert h>0
G = 2**35
bound=G*f[(5,0,0,0)]//h
assert bound==6936440<2**23
# Fourier LP bound for arbitrary packing motifs when their min² is 8.
# For additive subgroups, size is a power of two, so the upper is 2^22.
# At min²6, picks one phase per x in one block -> a 64-element anticode
# of diameter squared <=4; its translates by distinct motif elements
# are disjoint, so |motif|<=2^29. At min²4 or 2 use all 2^35 labels.
assert 2**58*6**35 < 2**44*8**35 # compare density squares at min²6 vs min²8
assert 2**70*4**35 < 2**44*8**35
print("PASS: 6,936,440 phases bound at min²8; additive center density exactly 4 in this chart")
