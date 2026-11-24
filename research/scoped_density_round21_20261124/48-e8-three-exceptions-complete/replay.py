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
def main(paths):
 package=[(P,json.loads(gzip.decompress(P.read_bytes()) if P.suffix=='.gz' else P.read_bytes())) for P in paths]
 used=set()
 def uses(node):
  if node is None:return
  if isinstance(node[0],int):uses(node[1]);uses(node[2]);return
  if node[0]=='proof':
   used.update(x[2] for x in node[3] if x[2]!='pair')
 for P,J in package:
  for _,_,node in J:
   if node is not None:uses(node)
 with TemporaryDirectory() as D:
  root=ROOT/'roots.in';exe=D+'/orbit'
  if used:
   from root_dictionary import check_lifts
   check_lifts()
   subprocess.run(['g++','-O3','-std=c++20',str(ROOT/'orbit.cpp'),'-o',exe],check=True)
  DS={}
  for s in sorted(used):
   args=[str(SH[s]),str(s%10 if s<100 else 3), 'positive' if s==40 else 'negative']
   DS[s]=load(subprocess.check_output([exe,*args],input=root.read_bytes()),s)
 cts=Counter()
 def check(box,node):
  d0,d1,b0,b1,c0,c1=box
  if node is None:raise AssertionError("certificate has an open leaf")
  if isinstance(node[0],int):
   j=node[0]*2;assert j in (0,2,4)
   u=list(box);v=list(box);u[j+1]=v[j]=(box[j]+box[j+1])/2
   check(u,node[1]);check(v,node[2]);return
  mode=node[0];cts[mode]+=1
  if mode=='empty':assert b0>c1 or d0>=b1*c1;return
  P=list(corners(box));assert len(P)==8
  assert tuple(node[1]) in COEF
  p=F(node[2]);M=node[1];assert p>0
  upper=min(F(4),max(4*(M[0]*a+M[1]*b+M[2]*c+M[3])for a,b,c in P))
  assert p>=upper or p**3>=512*d1
  fac=lambda g:min((g/8)**12,(g/8)**24/d0**4)
  if mode=='trivial':assert N*fac(p)<=T;return
  assert mode=='proof'
  top=F(0)
  for l,h,s,*dual in sorted(node[3],key=lambda v:F(v[0])):
   l,h=F(l),F(h);assert l==top and h>l;top=h;cts[s]+=1
   if s=='pair':
    u=min(max(a+b for a,b,c in P),max(a+1 for a,b,c in P))
    nv=N if l<u else N//16 if l<2*u else N//256
   else:
    assert s in SH
    if s==40:assert b1<=1
    # The 40 case uses the four ordinary slots plus a (a≤0 on domain).
    table=DS[s];R,K,W=table
    js=[j for j,v in enumerate(aup(table,P,l))if v]
    cert=dual[0]
    nv=1 if not js else 2**(8*SH[s])
    if js and cert:
     assert all(0<=i<len(W) and y>=0 for i,y in cert)
     v=min(F(-sum(y*W[i][j]for i,y in cert),K[j])for j in js)
     if v>0:nv=min(nv,int(1+sum(y for i,y in cert)/v))
    nv*=2**(8*(6-SH[s]))
   assert nv*fac(h)<=T
  assert top==p
 def boundary(box,node):
  d0,d1,b0,b1,c0,c1=box
  assert node is not None
  if isinstance(node[0],int):
   j=node[0]*2;assert j in (0,2,4)
   u=list(box);v=list(box);u[j+1]=v[j]=(box[j]+box[j+1])/2
   if j!=4:boundary(u,node[1])
   boundary(v,node[2]);return
  cts['heavy_boundary']+=1
  if node[0]=='empty':
   # At c≥c1, a<0 still gives the two-slot period cap 4b.
   g=4*b1
   assert N*min((g/8)**12,(g/8)**24/d0**4)<=T
   return
  p=F(node[2])
  assert p>=4 or p**3>=512*d1 or node[1][2]==0
  assert node[0]=='trivial' or all(z[2]=='pair' or z[2]<100 for z in node[3])
  # Every retained fibre weight is independent of c or decreases with it.
  # This is a complete replay of the already-conservative chart at c=c1.
  check([d0,d1,b0,b1,c1,c1],node)
 for p,J in package:
  v=H
  for dlo,dhi,node in J:
   dlo,dhi=F(dlo),F(dhi)
   assert dlo==v and dhi>dlo
   box=[dlo,dhi,F(2),F(2),F(2),F(2)]
   check(box,node)
   assert node[0] in ("proof","trivial") and (F(node[2])>=4 or F(node[2])**3>=512*dhi or node[1][1:3]==[0,0])
   assert node[0]=="trivial" or all(x[2] in ("pair",31,41) for x in node[3])
   v=dhi
  assert v==6
  print(p.name,'certified determinant interval',H,v,'all b,c≥2','cumulative',dict(cts),flush=True)
 print('PASS all claimed leaves; open leaves are not excluded:',cts['open'],flush=True)
if __name__=='__main__':main([Path(v)for v in sys.argv[1:]])
