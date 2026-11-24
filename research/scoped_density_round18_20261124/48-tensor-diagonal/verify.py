#!/usr/bin/env python3
import json
import subprocess
from pathlib import Path
from struct import unpack_from
from fractions import Fraction as F
from tempfile import TemporaryDirectory
from collections import Counter
ROOT=Path(__file__).resolve().parent
TARGET=F(38246144,2401)
HERM=F(389,2000)
B6=38246144
N0=2**48
def check_lifts():
 from math import prod
 def q(x):return (x&(x>>4)).bit_count()%2
 def b(x,y):return q(x^y)^q(x)^q(y)
 def dot(x,y):return sum(a*b for a,b in zip(x,y))
 R=[tuple(map(int,l.split())) for l in (ROOT/'roots.in').read_text().splitlines()]
 A=[x for x in range(256) if q(x)]
 assert len(R)==256 and len(A)==120 and all(len(r)==8 for r in R)
 for x,v in enumerate(R):
  if x not in A:
   assert not any(v)
   continue
  assert dot(v,v)==8 and len({k%2 for k in v})==1 and sum(v)%4==0
  for y in A:
   z=dot(v,R[y])
   assert z%4==0 and (z//4)%2==b(x,y)
 span={0}
 for x in A:
  span|={x^u for u in span}
 assert len(span)==256
 base=[1,16,2,32,4,64,8,128]
 order=[135,64,35,16,9,4,2,1]
 extra=[21,25,29,17]
 assert all(x in A for x in extra)
 def tr(x,r):return x^r if b(x,r) else x
 def comp(x):
  for r in extra:x=tr(x,r)
  return x
 assert all(comp(j)==j for j in base[:5])
 for i, x in enumerate(base):
  ts=[r for r in A if all(not b(r,t) for t in base[:i])]
  orbit={x}; stack=[x]
  for u in stack:
   w=[tr(u,r) for r in ts]
   if i==5:w.append(comp(u))
   for v in w:
    if v not in orbit:orbit.add(v);stack.append(v)
  assert len(orbit)==order[i]
 size=prod(order)
 assert size==prod(((1<<(2*m-1))+(1<<(m-1))-1)*(1<<(2*m-2)) for m in range(1,5))==348364800
 print('E8 root dictionary and surjectivity of liftable orthogonal group PASS',flush=True)

def load(raw,code):
    c,r=unpack_from('<II',raw); N=[]; L=[]; pos=8
    assert c==r
    for i in range(c):
        n,s=unpack_from('<qH',raw,pos); pos+=10
        assert n>0
        N.append(n)
        w=unpack_from('<%dH'%(3*s),raw,pos); pos+=6*s
        L.append([w[k:k+3] for k in range(0,3*s,3)])
    W=[unpack_from('<%dq'%c,raw,pos+8*c*i) for i in range(r)]
    assert pos+8*c*r==len(raw) and sum(N)+1==2**(8*(code//10))
    assert all(sum(row)==-1 for row in W)
    return N,L,W

def fibre_cap(N,W,js,k):
    assert isinstance(k,list) and all(isinstance(t,list) and len(t)==2 for t in k)
    assert all(type(i)==int and 0<=i<len(W) and type(y)==int and y>0 for i,y in k)
    assert len({i for i,_ in k})==len(k)
    if not js:
        return F(1)
    z=min(F(-sum(y*W[i][j] for i,y in k),N[j]) for j in js)
    assert z>0
    return 1+F(sum(y for _,y in k),z)

def uniform6(raw):
    N,*W=[tuple(map(int,line.split())) for line in raw.decode().splitlines()]
    assert len(N)==347 and len(W)==449 and all(n>0 for n in N)
    d=json.loads((ROOT/'old8.json').read_text())
    k=[[int(i),y] for i,y in d['numerators'].items()]
    z=fibre_cap(N,W,range(len(N)),k)
    assert z==F(d['cap']) and int(z)*256==B6
    print("Uniform six-slot fibre cap",int(z),flush=True)

def det(a,b):return a+b+5*a*b

def period(a,b):
    return min([F(8)]+[4*(a*(n+5*r)**2+b*n*n+5*r*r)
                for n in range(-8,9) for r in range(9) if (n,r)!=(0,0)])

BOXES={
    'neg':(F(-1,6),F(0),F(0),F(1)),
    'mat':(F(0),F(1),F(0),F(1)),
    'mneg':(F(-1,5),F(0),F(1),F(1)),
    'mpos':(F(0),F(1),F(1),F(1)),
    'high':(F(1),F(2),F(0),F(0)),
    'highinf':(F(2),F(2),F(0),F(0))}
SHAPES=(22,31,32,40,41,42)

def check_tree(chart,DS,tr,box,stats):
    ax,ay,bx,by=box
    if isinstance(tr,dict):
        assert len(tr)==1
        ((axis,(lt,rt)),)=tr.items()
        assert axis in ('a','b')
        if axis=='a':
            assert ax<ay
            m=(ax+ay)/2
            check_tree(chart,DS,lt,(ax,m,bx,by),stats)
            check_tree(chart,DS,rt,(m,ay,bx,by),stats)
        else:
            assert bx<by
            m=(bx+by)/2
            check_tree(chart,DS,lt,(ax,ay,bx,m),stats)
            check_tree(chart,DS,rt,(ax,ay,m,by),stats)
        return
    assert isinstance(tr,list)
    stats['tiles']+=1
    if chart in ('mat','neg'):
        if ax>by or det(ay,by)<=HERM:
            assert not tr
            return
        D0=det(ax,bx)
        p=min(F(7),period(ay,by))
    elif chart in ('mneg','mpos'):
        D0=det(ax,F(1))
        p=min(F(7),4*(1+ay),20*(1+5*ay))
    elif chart=='high':
        D0=det(ax,ax)
        p=min(F(8),7*ay)
    else:
        assert chart=='highinf'
        D0=det(F(2),F(2))
        p=F(8)
    D0=max(D0,HERM)
    if p<=0:
        assert not tr
        return
    inf=F(10**8)
    z={s:[min((h+l*ay+m*by for h,l,m in L),default=inf)
          for L in ds[1]] for s,ds in DS.items()}
    cuts=sorted({F(0),p}|{y for q in z.values() for y in q if 0<y<p})
    i=0
    for lo,hi in zip(cuts,cuts[1:]):
        f=(hi/8)**24/D0**4
        if N0*f<=TARGET: continue
        s,k=tr[i]; i+=1
        N,_,W=DS[s]
        nd=fibre_cap(N,W,[j for j,y in enumerate(z[s]) if y>lo],k)
        nd=min(int(nd),2**(8*(s//10)))
        assert nd*2**(8*(6-s//10))*f<=TARGET,(chart,box,lo,hi,nd)
        stats['dual bands']+=1
    assert i==len(tr)

def main():
    # Exact common caps for metrics with small determinant and the high-g windows.
    assert F(20,3)**12*HERM**8<TARGET
    assert F(B6,7**4)==TARGET<F(3,2)**24
    check_lifts()
    with TemporaryDirectory() as tmp:
        exe=tmp+'/orbits'
        subprocess.run(['g++','-std=c++20','-O3',str(ROOT/'orbit.cpp'),'-o',exe],check=True)
        roots=(ROOT/'roots.in').read_bytes()
        with open(Path(tmp)/'progress.log','wb') as log:
            def run(*args):
                return subprocess.check_output([exe,*args],input=roots,stderr=log)
            uniform6(run('uniform8'))
            ds={(neg,s):load(run(str(s//10),str(s%10), 'negative' if neg else 'positive'),s)
                for neg in (False,True) for s in SHAPES if not (s==40 and neg)}
        for chart,box in BOXES.items():
            S= SHAPES if chart in ('mat','neg') else ((40,) if chart=='highinf' else (31,40,41))
            m={s:ds[chart in ('neg','mneg') and s!=40,s] for s in S}
            stats=Counter()
            check_tree(chart,m,json.loads((ROOT/f'proof_{chart}.json').read_text()),box,stats)
            print(chart,dict(stats),'PASS',flush=True)
    print('Full arbitrary half-period motif cap:',TARGET)

if __name__=='__main__': main()
