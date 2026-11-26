"""Periodic equality in the canonical projected 22D Leech grid; no dependencies."""
from collections import Counter
from itertools import combinations
raw=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[2,-2,0,0,0,-2,2,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,2],[0,2,-2,0,0,0,2,2,-2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,2],[2,0,-2,0,0,0,2,0,-2,2,0,0,0,0,2,0,2,0,0,0,0,0,0,2],[-4,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[4,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,-1,1,-1,1,-1,3,1,-1,1,1,1,-1,1,-1,-1,1,-1,-1,1,-1,1,1,1],[1,1,-1,1,-1,-1,3,1,-1,1,1,1,-1,1,-1,1,1,1,-1,1,1,1,1,1],[1,-1,-1,-1,-1,-1,3,1,-1,-1,1,1,1,1,1,-1,1,1,-1,-1,-1,1,1,1],[-1,1,1,-1,-1,1,3,1,-1,1,1,1,1,1,-1,-1,1,1,1,-1,1,1,1,1],[1,-1,1,-1,-1,1,3,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,1,1,1],[-1,1,-1,1,1,1,3,1,-1,1,-1,1,1,1,1,-1,1,1,-1,1,-1,1,1,1]]
alternatives=[[0,1,[2,-2,0,0,0,2,2,2,0,2,0,-2,0,0,0,0,0,0,0,0,0,0,0,2]],[1,1,[2,-2,0,0,0,2,2,2,0,2,0,-2,0,0,0,0,0,0,0,0,0,0,0,2]],[0,1,[2,-2,0,2,2,0,0,2,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2]],[1,1,[2,-2,0,2,2,0,0,2,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2]],[0,1,[2,-2,0,-2,-2,0,0,2,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2]],[1,1,[2,-2,0,-2,-2,0,0,2,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2]],[0,1,[2,-2,0,0,0,0,0,0,-2,2,2,0,0,2,0,0,0,0,2,0,0,0,0,2]],[1,1,[2,-2,0,0,0,0,0,0,-2,2,2,0,0,2,0,0,0,0,2,0,0,0,0,2]],[0,1,[2,-2,0,0,0,0,0,0,-2,2,2,0,0,-2,0,0,0,0,-2,0,0,0,0,2]],[0,1,[2,-2,0,0,0,2,0,0,-2,0,0,0,2,0,2,0,0,0,0,2,0,0,0,2]],[0,1,[2,-2,0,0,0,-2,0,0,-2,0,0,0,2,0,2,0,0,0,0,-2,0,0,0,2]],[0,1,[2,-2,0,0,0,0,2,0,-2,0,2,2,0,0,0,0,2,0,0,2,0,0,0,0]],[1,1,[2,-2,0,0,0,0,2,0,-2,0,2,2,0,0,0,0,2,0,0,2,0,0,0,0]],[0,1,[2,-2,0,0,2,0,2,0,-2,2,0,0,2,0,0,0,0,0,0,0,2,0,0,0]],[0,1,[2,-2,0,0,-2,0,2,0,-2,2,0,0,2,0,0,0,0,0,0,0,-2,0,0,0]],[0,1,[2,-2,0,0,0,0,2,0,-2,0,0,0,0,0,0,2,0,2,0,0,0,2,0,2]],[0,1,[2,-2,0,0,0,0,2,0,0,0,0,0,2,-2,0,0,2,0,0,0,0,0,2,2]],[0,1,[3,-1,1,1,-1,-1,1,1,-1,1,1,-1,1,1,1,-1,1,1,-1,1,1,1,1,1]],[0,1,[1,-3,-1,-1,1,-1,1,1,-1,1,1,-1,1,1,1,1,1,-1,-1,1,-1,1,1,1]],[0,1,[1,-3,-1,1,-1,1,1,1,-1,1,1,1,1,1,1,-1,1,-1,-1,-1,1,-1,1,1]],[0,1,[3,-1,1,-1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,1,1]],[0,2,[0,2,-2,0,0,0,2,2,-2,0,-2,0,-2,0,0,0,0,0,0,0,0,0,0,2]],[0,2,[0,2,-2,0,0,0,0,0,0,2,-2,2,0,0,2,0,0,0,0,2,0,0,0,2]]]
G=set()
for m in range(4096):
    c=0
    for j in range(12):
        if (m>>j)&1: c ^= 0xae3<<j
    G.add(c|((c.bit_count()%2)<<23))
assert Counter(c.bit_count() for c in G)=={0:1,8:759,12:2576,16:759,24:1}
def leech(x):
    m=x[0]%2
    return (len(x)==24 and all((u-m)%2==0 for u in x)
      and sum((((u-m)//2)%2)<<i for i,u in enumerate(x)) in G and (sum(x)-4*m)%8==0)
def W(x):return (sum(x[:3]),*x[3:])
def d24(x,y):return (x[0]-y[0])**2 + 3*sum((a-b)**2 for a,b in zip(x[1:],y[1:]))
assert all(map(leech,raw))
w=list(map(W,raw))
def profile(v):return Counter(d24(x,y) for x,y in combinations(v,2))
assert len(set(w))==12 and profile(w[:4])=={72:6} and max(profile(w))==88
B=[];sg=8559;assert sg in G and sg%8 in (0,7)
def flip(x,flag):return [t*(-1 if flag and (sg>>j)&1 else 1) for j,t in enumerate(x)]
for flag,k,alt in alternatives:
    base=[flip(x,flag) for x in raw[:4]]
    rep=flip(alt,flag);new=base[:];new[k]=rep
    assert all(map(leech,base+new))
    assert profile(list(map(W,base)))=={72:6} and profile(list(map(W,new)))=={72:6}
    assert set(map(W,base))&set(map(W,new))==set(W(base[i]) for i in range(4) if i!=k)
    diff=[a-b for a,b in zip(rep,base[k])]
    assert leech(diff) and diff[0]==diff[1]==diff[2]
    B.append(diff[2:])
# Integer row operations (swaps, minus integer multiples of rows, sign) preserve the row lattice.
def index(A):
    A=[r[:] for r in A];rank=0;det=1
    for j in range(22):
        pivot=next((k for k in range(rank,len(A)) if A[k][j]),None)
        if pivot is None:continue
        A[rank],A[pivot]=A[pivot],A[rank]
        for k in range(rank+1,len(A)):
            while A[k][j]:
                q=A[rank][j]//A[k][j]
                A[rank],A[k]=A[k],[x-q*y for x,y in zip(A[rank],A[k])]
        det*=abs(A[rank][j]);rank+=1
    assert rank==22
    return det
assert len(B)==23 and index(B)==2**34
def label(x):
    p=(x[0]-x[1])//2;q=(x[1]-x[2])//2
    return ((4*p+2*q)%12,(2*p+4*q)%12)
def add(a,b,s=1):return tuple((x+s*y)%12 for x,y in zip(a,b))
U={label(x): d24(W(x),W(raw[0])) for x in raw}
assert len(U)==12 and Counter(U.values())=={0:1,64:2,72:3,88:6}
allowed={frozenset(C) for C in combinations(U,3) if all(U[add(a,b,-1)]==88 for a,b in combinations(C,2))}
ref=[(0,0),label([2,0,0]),label([0,-2,0])]
orbit={frozenset(add(a,b,s) for b in ref) for a in U for s in [1,-1]}
assert allowed==orbit and len(orbit)==24
print('tetrahedron / 12-point profiles:',profile(w[:4]),profile(w))
print('forced sublattice raw determinant:',index(B),'rank 22, admissible color triples:',len(allowed))
