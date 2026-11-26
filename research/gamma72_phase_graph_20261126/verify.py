"""Exact audit, including a complete Leech shell reconstruction from the existing certificate."""
from pathlib import Path
from itertools import combinations
import json, runpy
import numpy as np
P=Path(__file__).resolve().parent
d=runpy.run_path(str(P.parent/'gamma72_depth6_20261122'/'verify.py'))
F=d['F']; A=d['A']; B=d['B']
def pack(r):
    return sum((int(x)&1)<<i for i,x in enumerate(r))
basis={}
for row in A:
    v=pack(row)
    while v:
        j=v.bit_length()-1
        if j not in basis: basis[j]=v; break
        v^=basis[j]
piv=sorted(basis,reverse=True)
assert len(piv)==12
def coordA(v):
    out=0
    for i,p in enumerate(piv):
        if v>>p&1: v^=basis[p];out|=1<<i
    assert not v
    return out
def coordB(v):
    w=pack(F@np.array([v>>i&1 for i in range(24)],dtype=np.int64))
    return sum(((basis[p]&w).bit_count()%2)<<i for i,p in enumerate(piv))
# Thus the polarized pairing compressed to 12 coordinates really is dot product.
M=[set() for _ in range(4096)]
for a,row in d['M4'].items():
    M[coordA(a)]={coordB(b) for b in row}
assert not M[0] and all(len(s)==24 for s in M[1:])
av=np.arange(4096,dtype=np.uint16)
pc=np.fromiter((v.bit_count()&1 for v in range(4096)),dtype=np.uint8)
ann=[set() for _ in range(4096)]
for a in range(1,4096):
    z=min(M[a]);ok=np.ones(4096,dtype=bool)
    for x in M[a]:
        ok &= pc[(x^z)&av]==0
    ann[a]=set(int(y) for y in av[ok])
assert {n:sum(len(v)==n for v in ann) for n in (4,16,64)}=={4:3380,16:650,64:65}
R=set();before=0
for a in range(1,4096):
    for b in ann[a]:
        c=a^b
        if not a<b<c or a not in ann[b] or a not in ann[c]:continue
        before+=1
        x,y,z=(min(M[w]) for w in (a,b,c))
        u=(a&y).bit_count()%2;v=(b&x).bit_count()%2;w=(c&x).bit_count()%2
        if (u,v,w,w)!=((a&z).bit_count()%2,(b&z).bit_count()%2,(c&y).bit_count()%2,u^v):continue
        # Keep the SUPerset of all possible deeply separated projected colors; ignoring 4+4 cannot generate a false no-clique conclusion.
        R.update((r<<12)|s for r in (a,b,c) for s in (a,b,c) if r!=s)
assert before==3315 and len(R)==5850
rv=np.array(sorted(R),dtype=np.uint32); flag=np.zeros(1<<24,dtype=bool);flag[rv]=True
for x in rv:
    N=rv[flag[rv^x]]
    assert len(N)==14
    for y in N:
        assert np.array_equal(N[flag[N^y]],np.array([x^y],dtype=np.uint32))
print('PASS all 5850 projected candidates; 14 neighbors each, link exactly seven disjoint edges')
def dep(a,b,total):
    if not (a or b or total):return 0
    if not (a and b and a!=b):
        u=a or b
        if total in M[u]:return 2
        if (u&total).bit_count()%2:return 3
        return 4
    tri=(a,b,a^b)
    for j,k in combinations(tri,2):
        if any((total^x) in M[k] for x in M[j]):return 4
    for j in tri:
        for k in tri:
            if j!=k and any(((total^x)&k).bit_count()%2 for x in M[j]):return 5
    return 6
def tab(raw):
    got=[];bound=[];Acols=[]
    for i in range(1<<len(raw)):
        v=np.array([i>>j&1 for j in range(len(raw))],dtype=np.int64)@raw%2
        u,w,s=d['syndrome'](v);a,b,t=coordA(u),coordA(w),coordB(s)
        got.append(dep(a,b,t));Acols.append((a,b))
        if not i:bound.append(0)
        elif not(a and b and a!=b):bound.append(4)
        else:bound.append(6 if ((a<<12)|b) in R else 5)
    return got,bound,Acols
payload=json.loads((P/'phases.json').read_text())
other=json.loads((P.parent/'gamma72_root_tilt_20261122'/'phases.json').read_text())
raw=np.array(payload['rows'],dtype=np.int64)
assert raw.shape==(7,72) and np.array_equal(raw,raw%2)
x,upper,asnew=tab(raw);old,_,asold=tab(np.array(other['old_D4']+other['new'],dtype=np.int64))
assert x==payload['depth_by_mask']==upper and asold==asnew
for k,count in ((5,2),(6,6),(7,24)):
    dif=[x[i]-old[i] for i in range(1<<k)]
    assert set(dif)=={0,2} and dif.count(2)==count
print('PASS all 128 parity depths simultaneously attain the maximum compatible with their A-projections; improved by 2 on 2,6,24 masks in ranks 5,6,7')
