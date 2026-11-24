"""Replay rational harmonic-shell duals, no numerical solver."""
from itertools import product, combinations_with_replacement
from fractions import Fraction as F
from pathlib import Path
import json

def exps(s):
    return [x for x in product(range(s+1),repeat=3) if sum(x)==s]

def lap(C,K):
    ans={}
    for e,a in C.items():
        for i in range(3):
            for j in range(3):
                f=list(e)
                b=f[i];f[i]-=1
                b*=f[j];f[j]-=1
                if b>0 and min(f)>=0:
                    g=tuple(f)
                    ans[g]=ans.get(g,0)+a*b*K[i][j]
    return ans

def harmonic(e,K,d):
    C={e:1};h=C.copy();f=F(1)
    for j in range(1,d//2+1):
        C=lap(C,K)
        f*=-F(8,2*j*(88+2*d-2*j-2))
        h.update({a:b*f for a,b in C.items()})
    return h

def val(C,x):
    return sum(a*x[0]**e[0]*x[1]**e[1]*x[2]**e[2] for e,a in C.items())

def add(C,H,u):
    for e,v in H.items():
        C[e]=C.get(e,0)+v*u

dct=json.loads(Path(__file__).with_name('duals.json').read_text())
known=set()
A=168498000
lo=A
for a,b,c in combinations_with_replacement(range(5),3):
    for u,v,w in ([(-a,-b,-c),(a,b,c)] if a else [(a,b,c)]):
        if u+v+w < -8:
            continue
        D=512-8*(u*u+v*v+w*w)+2*u*v*w
        if D>384:
            continue
        assert D>0
        k=f'{u},{v},{w}'
        known.add(k)
        g=[[8,u,v],[u,8,w],[v,w,8]]
        adj=[[64-w*w,v*w-8*u,u*w-8*v],
             [v*w-8*u,64-v*v,u*v-8*w],
             [u*w-8*v,u*v-8*w,64-u*u]]
        assert [[sum(x*y for x,y in zip(row,col)) for col in zip(*adj)] for row in g] == [[D*int(i==j) for j in range(3)] for i in range(3)]
        dual=list(map(F,dct[k]))
        H=[harmonic(e,g,d) for d in [2,6] for e in exps(d)]
        assert len(dual)==len(H)+1
        P={(0,0,0):dual[0]}
        for z,h in zip(dual[1:],H):
            add(P,h,z)
        for v0 in product(range(-4,5),repeat=3):
            if sum(x*y*a for x,row in zip(v0,adj) for y,a in zip(v0,row))<=8*D:
                assert val(P,v0)<=int(v0==(0,0,0)),(k,v0)
        lower=A*dual[0]-2*sum(val(P,z) for z in g)
        assert lower>1000000
        lo=min(lo,lower)
assert known==set(dct)
assert lo>1033558
print('PASS exact positive zero-label shell for all',len(known),'Gram shapes')
print('universal count lower bound',lo.numerator//lo.denominator)
