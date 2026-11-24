#!/usr/bin/env python3
import gzip,json,math,subprocess
from collections import Counter
from functools import cache
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
import base
ROOT=Path(__file__).resolve().parent
TARGET=F(16800)
H=F(11,4000)
N0=2**48
shapes=(22,31,32,40,41,42)
def det(a,b,c):return a*b+a*c+b*c+4*a*b*c
period=sorted({((u+v+4*t)**2,u*u,v*v,4*t*t)
 for u,v,t in product(range(-5,6),repeat=3) if (u,v,t)!=(0,0,0)})
def per(a,b,c):
 den=math.lcm(a.denominator,b.denominator,c.denominator)
 x,y,z=int(a*den),int(b*den),int(c*den)
 return 4*F(min(h*x+j*y+k*z+l*den for h,j,k,l in period),den)
@cache
def affine(neg,s,a,b):
 # Every stored lift length is a homogeneous sum of slot weights.
 # Treating the third slot as either fixed or ordinary upper weight 1
 # can only increase these witness lengths (when its value is <=1).
 L=DS[neg and s!=40,s][1]
 return tuple(min((h+l*a+m*b for h,l,m in vals),default=F(10**8)) for vals in L)
@cache
def cap(neg,s,js,dual):
 N,_,W=DS[neg and s!=40,s]
 return min(int(base.fibre_cap(N,W,js,[list(x) for x in dual])),2**(8*(s//10)))
def bands(neg,shape,a,b,p,D,tr,counts):
 if p<=0:
  assert tr==[]
  return
 vals={s:affine(neg,s,a,b) for s in shape}
 cuts=sorted({F(0),p}|{v for z in vals.values() for v in z if 0<v<p})
 i=0
 for lo,hi in zip(cuts,cuts[1:]):
  fac=(hi/8)**24/D**4
  if N0*fac<=TARGET:continue
  s,k=tr[i];i+=1
  assert s in shape
  z=vals[s];U=tuple(j for j,y in enumerate(z) if y>lo)
  v=cap(neg,s,U,tuple(map(tuple,k)))
  assert v*2**(8*(6-s//10))*fac<=TARGET,(s,hi,D,v)
  counts['dual bands']+=1
 assert i==len(tr)
def check(chart,box,tr,ct):
 if isinstance(tr,dict):
  assert len(tr)==1
  ((key,(lt,rt)),)=tr.items()
  i=2*('abc'.index(key) if chart!='neg3' else int(key))
  assert 0<=i<len(box) and box[i]<box[i+1]
  m=(box[i]+box[i+1])/2
  a,b=list(box),list(box)
  a[i+1]=m;b[i]=m
  check(chart,a,lt,ct);check(chart,b,rt,ct)
  return
 assert isinstance(tr,list)
 ct['tiles']+=1
 if chart=='mat3':
  a,A,b,B=box
  if a>B or A+B<=F(3,5):
   assert tr==[]
   return
  D=max(det(a,b,b),F(27,100))
  bands(False,shapes,A,B,min(F(7),4*(A+B)),D,tr,ct)
 elif chart=='neg3':
  a,A,b,B,c,C=box
  if b>C:
   assert tr==[]
   return
  d=[det(*x) for x in product([a,A],[b,B],[c,C])]
  if max(d)<=H:
   assert tr==[]
   return
  D=max(H,min(d))
  bands(True,shapes,A,B,min(F(7),per(A,B,C)),D,tr,ct)
 else:
  assert chart=='high3'
  a,A,_,_=box
  D=det(a,a,a)
  bands(False,(30,31,41),A,F(0),min(F(8),7*A),D,tr,ct)
def run():
 global DS
 # Hermite cusp and near-origin shortcuts.
 assert N0*H**4<TARGET
 assert 2**24*F(4,3)**4*F(3,5)**16<TARGET
 assert TARGET<F(3,2)**24
 base.check_lifts()
 with TemporaryDirectory() as tmp:
  exe=tmp+'/orbit'
  subprocess.run(['g++','-O3','-std=c++20',str(ROOT/'orbit.cpp'),'-o',exe],check=True)
  roots=(ROOT/'roots.in').read_bytes()
  def regenerate(*args):return subprocess.check_output([exe,*args],input=roots,stderr=subprocess.DEVNULL)
  base.uniform6(regenerate('uniform8'))
  DS={(neg,s):base.load(regenerate(str(s//10),str(s%10),'negative' if neg else 'positive'),s)
      for neg in (False,True) for s in shapes+(30,) if not(neg and s in (30,40))}
 # For heavy c with a<=1, transfer the prior two-exception windows g<=7.
 assert base.TARGET<=TARGET
 assert 2**24*base.HERM**20<TARGET
 for name,which in [('mat',shapes),('mpos',(31,40,41))]:
  stats=Counter()
  tab={s:DS[False,s] for s in which}
  base.check_tree(name,tab,json.loads((ROOT/f'prior_{name}.json').read_text()),base.BOXES[name],stats)
  print('inherited',name,dict(stats),'PASS',flush=True)
 # Small Fourier caps for a>=2, or g>7*max(1,a) with at least two heavy slots.
 small=json.loads((ROOT/'small.json').read_text())
 for s,t,M in [(30,5,1606),(40,7,1243)]:
  nd=cap(False,s,tuple(j for j,d in enumerate(affine(False,s,F(0),F(0))) if d>t),tuple(map(tuple,small[str(s)])))
  assert nd<=M
 assert 1606*256**3/F(44)**4<TARGET and N0*F(5,8)**24/F(44)**4<TARGET
 assert 1243*256**2*F(49,432)**4<TARGET
 assert F(38246144,2401)<TARGET
 # All-light continuum and positive-heavy a in [1,2].
 boxes={
 'mat3':(F(0),F(1),F(0),F(1)),
 'neg3':(F(-1,6),F(0),F(0),F(1),F(0),F(1)),
 'high3':(F(1),F(2),F(0),F(0))}
 for name,box in boxes.items():
  stats=Counter()
  tr=json.loads(gzip.decompress((ROOT/f'{name}.json.gz').read_bytes()))
  check(name,box,tr,stats)
  print(name,dict(stats),'PASS',flush=True)
 print('Uniform 48D scoped center density cap:',TARGET)
if __name__=='__main__':run()
