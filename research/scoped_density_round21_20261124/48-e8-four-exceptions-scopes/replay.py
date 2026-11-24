#!/usr/bin/env python3
import json, gzip, math, sys, subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
from itertools import product
from collections import Counter
from struct import unpack_from
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parent
H=F(11,4000);N=2**48;T=16800
assert N*H**4<T and F(2**24,6**4)<T
SH={22:2,31:3,32:3,40:4,41:4,42:4,102:2,103:3,104:4}
ZERO=2**60
def load(raw,s):
 c,r=unpack_from('<II',raw);p=8;R=[];sz=[];n=3 if s<100 else 4
 for _ in range(c):
  z,k=unpack_from('<qH',raw,p);p+=10;assert z>0
  a=unpack_from('<%dH'%(n*k),raw,p);p+=n*k*2;sz.append(z)
  R.append(tuple(tuple(a[j:j+n])for j in range(0,len(a),n)))
 W=[unpack_from('<%dq'%c,raw,p+i*c*8)for i in range(r)]
 assert p+8*c*r==len(raw) and sum(sz)+1==2**(8*SH[s]) and all(sum(w)==-1 for w in W)
 return R,sz,W
COEF={((u+v+t)**2,u*u,v*v,(t//4)**2*4+(2*(t//4)+1)*(t%4))for u,v,t in product(range(-10,11),range(-10,11),range(-28,29)) if (u,v,t)!=(0,0,0)}
def corners(box):
 for D,b,c in product(box[:2],box[2:4],box[4:]):
  yield (D-b*c)/(b+c+4*b*c),b,c
def aup(table,points,t):
 tab,_,_=table
 HZ=10**8
 A=[(HZ,math.ceil(a*HZ),math.ceil(b*HZ),math.ceil(c*HZ))for a,b,c in points]
 return tuple(min((max(sum(u*v for u,v in zip(w,p))for p in A)for w in src),default=ZERO)>t*HZ for src in tab)
def main():
 from root_dictionary import check_lifts
 check_lifts()
 with TemporaryDirectory() as tmp:
  exe=tmp+'/orbit'
  subprocess.run(['g++','-O3','-std=c++20',str(ROOT/'orbit.cpp'),'-o',exe],check=True)
  tab=load(subprocess.check_output([exe,'3','1','negative'], input=(ROOT/'roots.in').read_bytes()),31)
 R,K,W=tab
 N=2**48
 for r in (8,24): assert max(F(1,4)**r,F(1,16)*F(1,2)**r,F(1,256))==F(1,256)
 assert (N//256)*F(3,4)**6*F(1,3500)**2<16800
 assert F(2**24,6**4)<16800
 D=F(128,1875)
 counts=Counter()
 for d0,d1,node in json.loads(gzip.decompress((ROOT/'proof.json.gz').read_bytes())):
  d0,d1=F(d0),F(d1)
  assert d0==D and d1>d0
  D=d1
  a=(d1-8)/36
  p,rows=node;p=F(p)
  assert p>0 and (p>=min(4,4*(1+a),12*(1+3*a)) or p**4>=3072*d1)
  fac=lambda g:min((g/8)**24/d0**4,F(3,4)**4*(g/8)**8)
  if not rows:
   counts['trivial']+=1
   assert N*fac(p)<=16800
   continue
  top=F(0)
  for u,v,s,*dual in rows:
   u,v=F(u),F(v)
   assert u==top and u<v
   top=v
   counts[s]+=1
   if s=='pair':
    cap=N if u<a+1 else N//16 if u<2*(a+1) else N//256
   else:
    assert s==31
    js=[j for j,x in enumerate(aup(tab,[(a,0,0)],u)) if x]
    cert=dual[0]
    c=1 if not js else 2**24
    if js and cert:
     assert all(0<=i<len(W) and y>=0 for i,y in cert)
     z=min(F(-sum(y*W[i][j]for i,y in cert),K[j])for j in js)
     if z>0:c=min(c,int(1+sum(y for i,y in cert)/z))
    cap=c*2**24
   assert cap*fac(v)<=16800
  assert top==p
 assert D==6
 lo=F(1,3500)
 while lo < F(128,1875):
  hi=min(F(128,1875),2*lo)
  P=[]
  for b,c,d in product([F(2),F(5,2)],repeat=3):
   a=(hi-b*c*d)/(b*c+b*d+c*d+3*b*c*d)
   P.append(4*(81*a+b+c+d+12))
  g=max(P)
  assert g>0 and (N//256)*min((g/8)**24/lo**4,F(3,4)**4*(g/8)**8)<=16800
  lo=hi
 print('PASS: D<=1/3500 or D>=128/1875, and all D on [2,5/2]^3',dict(counts))
if __name__=='__main__': main()
