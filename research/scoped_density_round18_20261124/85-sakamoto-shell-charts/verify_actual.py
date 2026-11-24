"""Exact graph check using only literal parent vectors in one 2A3 chart."""
from sys import argv
from itertools import product
from pathlib import Path
import gzip, json, hashlib

def dot(x,y): return sum(a*b for a,b in zip(x,y))
def key(b): return tuple((4*a-sum(b))%16 for a in b)
def diff(x,y): return tuple((a-b)%16 for a,b in zip(x,y))
zero=(0,0,0)
def load(prefix):
    path=Path(__file__).parents[0]
    if (path/'actual_input.json.gz').exists():
        data=json.load(gzip.open(path/'actual_input.json.gz','rt'))
        return data[prefix]
    return json.loads((path/(prefix+'_lattice.json')).read_text())['G'], json.loads((path/(prefix+'_chart_witness.json')).read_text())

def omega(A):
    best=0
    def col(P):
        order=[]
        bounds=[]
        d=0
        while P:
            d+=1
            Q=P
            while Q:
                b=Q&-Q
                u=b.bit_length()-1
                P^=b; Q^=b; Q&=~A[u]
                order.append(u); bounds.append(d)
        return order,bounds
    def rec(P,s):
        nonlocal best
        if not P:
            best=max(best,s); return
        order,ub=col(P)
        for u,k in zip(reversed(order),reversed(ub)):
            if s+k<=best:return
            rec(P&A[u],s+1)
            P^=1<<u
    rec((1<<len(A))-1,0)
    return best

for tag in (argv[1:] or ['L','L0']):
    G,D=load(tag)
    sha={'L':'f95d33297e2713f696b13d5e123d24b86f2dd3212ee514fea6cafc2548f11c06',
          'L0':'53cf44a2d5ab1af9b4c2d93ad8b54d84041f1a4871874e6d8be4ffcce2bb247a'}
    assert hashlib.sha256(json.dumps(G,separators=(',',':')).encode()).hexdigest()==sha[tag]

    def gmul(x): return [dot(x,c) for c in zip(*G)]
    C=D['normals']
    B=[gmul(v) for v in C]
    assert [[dot(v,w) for w in B] for v in C]==[[8 if i==j else 4 for j in range(3)] for i in range(3)]
    S={key(b) for b in product(range(16),repeat=3)}
    assert len(S)==256
    block={}
    for cls,h,w in D['witnesses']:
        b=[dot(w,v) for v in B]
        z=key(b)
        assert z==tuple(cls) and z in S
        q=dot(w,gmul(w))
        assert q in (8,10)
        q16=4*sum(a*a for a in b)-sum(b)**2
        assert h==128-16*q+q16
        if z!=zero:
            block[z]=max(block.get(z,0),h)
            nz=diff(zero,z)
            block[nz]=max(block.get(nz,0),h)
    limits = []
    for h in sorted({0,*block.values()}):
        if h>=128:continue
        if 256**2*(128-h)**85<=64*128**85:break
        T=sorted(z for z in S if z!=zero and block.get(z,0)<=h)
        A=[sum(1<<i for i,w in enumerate(T) if z!=w and block.get(diff(z,w),0)<=h) for z in T]
        m=1+omega(A)
        ok=m*m*(128-h)**85<=64*128**85
        limits.append((h,m,ok))
    print(tag,'literal blocked',sum(h>0 for h in block.values()),'limits',limits)
    assert all(m[2] for m in limits)
    print('CLOSED')
