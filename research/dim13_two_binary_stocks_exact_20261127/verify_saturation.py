"""Exact covering upper bound for the literal Julin stock."""
from fractions import Fraction as F
from itertools import combinations

N = 4096
ALL = N-1
rows = [0] + [sum(1 << (1+(i+j)%11) for j in (0,1,3,4,5,9))
              for i in range(11)]
hexads = {a^b for a,b in combinations(rows,2)}
hexads |= {ALL^x for x in hexads}
loose = {3<<(2*i) for i in range(6)}
J = hexads | loose | {ALL^x for x in loose}
assert len(J)==144 and min((a^b).bit_count() for a,b in combinations(J,2))==4

# After keeping k coordinates, cover the retained binary cube within r.
for k,r in [(5,0),(9,1),(11,2),(12,3)]:
    for I in combinations(range(12),k):
        mask = sum(1<<i for i in I)
        P = {j&mask for j in J}
        short = [sum(1<<i for i in p) for d in range(r+1)
                 for p in combinations(I,d)]
        covered = {x^y for x in P for y in short}
        assert len(covered) == 1<<k


def solve(A,b):
    m=len(A)
    M=[list(map(F,row))+[F(b[i])] for i,row in enumerate(A)]
    for i in range(m):
        j=next((j for j in range(i,m) if M[j][i]),None)
        if j is None:
            return
        M[i],M[j]=M[j],M[i]
        p=M[i][i]
        M[i]=[x/p for x in M[i]]
        for j in range(m):
            if j!=i:
                p=M[j][i]
                M[j]=[x-p*y for x,y in zip(M[j],M[i])]
    return [row[-1] for row in M]


H=[lambda a:sum(a[:7]), lambda a:sum(a[:3])+a[11],
   lambda a:a[0]+a[10]+a[11], lambda a:sum(a[9:])]
V=[tuple(int(i>=12-k) for i in range(12)) for k in range(13)]
value=[[f(v) for f in H] for v in V]
Best=F(0)
for s in range(1,5):
    for I in combinations(range(13),s):
        for K in combinations(range(4),s):
            M=[[1]*s]+[[value[i][K[0]]-value[i][j] for i in I]
                                           for j in K[1:]]
            x=solve(M,[1]+[0]*(s-1))
            if x is None or min(x)<0:
                continue
            a=[sum(u*V[i][j] for u,i in zip(x,I)) for j in range(12)]
            bound=sum((1-t)**2 for t in a)/4 + min(f(a) for f in H)
            Best=max(Best,bound)
assert Best == F(103,32)
assert Best+F(1,4)==F(111,32)<4
print("PASS: Julin projection covers and all-real 13D saturation")

