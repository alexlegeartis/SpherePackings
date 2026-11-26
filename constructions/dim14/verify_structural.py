"""Exact finite inputs for the 13D/14D rigidity proofs and 14D/15D stack comparisons."""
from pathlib import Path
from itertools import combinations, product
from runpy import run_path
d=run_path(str(Path(__file__).with_name('verify.py')))
C=sorted(d['C']); cosets=d['cosets']; H=cosets[0]
def syndrome(w):
    s=0
    for i in range(12):
        if (w>>i)&1: s^=16+i
    return s
assert {syndrome(c) for c in C}==set(range(16)) and syndrome(4095)==0
pair_supports=[(1<<i)|(1<<j) for i,j in combinations(range(12),2)]
support_sets=[]
for J in (C,H):
    assert all(sum((c>>i)&1 for c in J)*2==len(J) for i in range(12))
    A={x:{x^y for y in J if (x^y).bit_count()==4} for x in J}
    S=set().union(*A.values()); support_sets.append(S)
    assert all(any(w&p==p for w in A[x]) for x in J for p in pair_supports)
    for i in range(12):
        V={J[0]};Q=[J[0]]
        for x in Q:
            for v in A[x]:
                y=x^v
                if (v>>i)&1 and y not in V:
                    V.add(y); Q.append(y)
        assert len(V)==len(J)
assert len(support_sets[0])==495 and len(support_sets[1])==39
assert all(any((s&t)==1<<i for s in support_sets[0] for t in support_sets[1]) for i in range(12))
for i,j in combinations(range(12),2):
    for a,b in product(range(2),repeat=2):
        assert sum((x>>i)&1==a and (x>>j)&1==b for x in H)==32
depth=[min(x.bit_count() for x in C) for C in cosets]
assert depth==[0]+[2]*15+[1]*12+[3]*4
def fcc(x):
    return x[0]%2,(28*(x[1]%2))^((sum(i//2 for i in x))%2)
for x in product(range(4),repeat=3):
    if sum(x)%2: continue
    phase,p=fcc(x)
    for dx in product(range(-2,3),repeat=3):
        n=sum(i*i for i in dx)
        if sum(dx)%2 or not 0<n<8: continue
        other,q=fcc(tuple(a+b for a,b in zip(x,dx)))
        h=3 if phase!=other else depth[p^q]
        assert n+2*h>=8
print('PASS: exact Julin/Hamming contact combinatorics, difference covering inputs, and FCC tie.')
