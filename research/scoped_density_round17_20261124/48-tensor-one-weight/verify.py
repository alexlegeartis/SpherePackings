#!/usr/bin/env python3
import json
import sys
import subprocess
from pathlib import Path
from struct import unpack
from fractions import Fraction as F
from tempfile import TemporaryDirectory
from collections import Counter
ROOT=Path(__file__).resolve().parent
TARGET=F(3,2)**24
PER=(4,4); PER2=(24,144)
OLD=[0]+[2**48]*7+[38246144]
def load(blob):
 c,r=unpack('<II',blob[:8]); assert len(blob)==8+24*c+8*c*r
 cols=[unpack('<q8H',blob[8+24*i:32+24*i]) for i in range(c)]
 rows=[unpack('<%dq'%c,blob[8+24*c+8*c*i:8+24*c+8*c*(i+1)]) for i in range(r)]
 for z in rows: assert sum(z)==-1
 return cols,rows
def val(t,a):
 a=F(a)
 return t[0]+t[1]*a
def supports(cols):
 return [[(j,l) for j,l in enumerate(c[1:]) if l!=65535] for c in cols]
def cells(datas,negative):
 A0,A1=(F(-1,6),F(0)) if negative else (F(0),F(1))
 lines={(j,0) for j in range(9)}|(set((PER,PER2)) if negative else {PER})
 for cols,_ in datas.values():
  for ps in supports(cols): lines.update(ps)
 breaks={A0,A1}
 for h,l in lines:
  for j,k in lines:
   if l!=k and A0<F(j-h,l-k)<A1: breaks.add(F(j-h,l-k))
 bs=sorted(breaks); out=[]
 for a,b in zip(bs,bs[1:]):
  mid=(a+b)/2
  period=min((PER,PER2),key=lambda x:val(x,mid)) if negative else PER
  seq=sorted((x for x in lines if 0<=val(x,mid)<=val(period,mid)), key=lambda x:val(x,mid))
  assert seq[0]==(0,0) and seq[-1]==period
  for lo,hi in zip(seq,seq[1:]):
   assert val(lo,mid)<val(hi,mid)
   # At equality to lower endpoint use the preceding band. At either
   # endpoint of the a interval, the same excluded columns remain <=lo.
   gs=val(lo,mid)
   oldn=2**48 if negative else OLD[gs.numerator//gs.denominator+1]
   sets=tuple(tuple(j for j,p in enumerate(supports(cols))
       if all(val(x,mid)>gs for x in p)) for cols,_ in datas.values())
   fac=max((val(hi,x)/8)**24/(1+6*x)**4 if x>F(-1,6) else F(0) for x in (a,b))
   out.append((a,b,lo,hi,oldn,sets,fac))
 return out
def cap(cols,rows,js,k):
 assert all(0<=i<len(rows) and type(y)==int and y>0 for i,y in k)
 assert len(set(i for i,_ in k))==len(k)
 if not js:return 1
 z=min(F(-sum(y*rows[i][j] for i,y in k),cols[j][0]) for j in js)
 assert z>0
 return 1+F(sum(y for i,y in k),z)
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

def main():
 check_lifts()
 with TemporaryDirectory() as tmp:
  exe=tmp+'/weights'
  subprocess.run(['g++','-std=c++20','-O3',str(ROOT/'low_dim.cpp'),'-o',exe],check=True)
  ds={3:load(subprocess.check_output([exe,'3'],input=(ROOT/'roots.in').read_bytes()))}
  raw=subprocess.check_output([exe,'uniform8'],input=(ROOT/'roots.in').read_bytes()).decode()
  N,*W=[tuple(map(int,line.split())) for line in raw.splitlines()]
  assert len(N)==347 and len(W)==449
  d=json.loads((ROOT/'old8.json').read_text())
  k=[(int(i),y) for i,y in d['numerators'].items()]
  z=cap([(n,) for n in N],W,range(len(N)),k)
  assert z==F(d['cap']) and (z.numerator//z.denominator)*256==OLD[8]
  print('Uniform g>7 fixed-exceptional-label cap',OLD[8],flush=True)
  neg={h:load(subprocess.check_output([exe,str(h),'negative'],input=(ROOT/'roots.in').read_bytes())) for h in (2,3)}
 check_cells(ds,False)
 check_cells(neg,True)
def check_cells(ds,negative):
 dd=json.loads((ROOT/('negative.json' if negative else 'positive.json')).read_text())
 H=list(ds)
 jobs=cells(ds,negative); used=Counter(); remain=[]; maximum=F(0)
 for a,b,lo,hi,oldn,sets,fac in jobs:
  best=oldn*fac
  hid='uniform'
  for h,js in zip(H,sets):
   key=str(h)+':'+';'.join(map(str,js))
   if key not in dd:continue
   n=cap(*ds[h],js,dd[key]); n=(n.numerator//n.denominator if isinstance(n,F) else n)
   bound=n*(1<<(8*(6-h)))*fac
   if bound<best: best=bound; hid=str(h+1)
  used[hid]+=1
  if best>maximum: maximum=best; max_window=(a,b,lo,hi,hid)
  if best>=TARGET:remain.append((str(a),str(b),lo,hi,float(best)))
 print('Intervals, windows, witnesses',len(set((a,b) for a,b,*_ in jobs)),len(jobs),len(dd))
 print('Winning cap sources',dict(used),'maximum',float(maximum))
 print('Unresolved',len(remain))
 for v in remain:print(v)
 assert not remain
 if negative:
  assert len(jobs)==701 and len(dd)==8 and maximum==F(1213025622610333925376,79792266297612001)
 else:
  assert len(jobs)==96 and len(dd)==18 and maximum==F(38246144,2401)
 assert maximum<TARGET
 print('Full negative bound' if negative else 'Full nonnegative bound',maximum,float(maximum),max_window)
if __name__=='__main__': main()
