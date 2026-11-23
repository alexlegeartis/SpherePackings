"""Exact Selling-cut enumeration and 39D density comparisons."""
import numpy as np
from itertools import combinations,product
from fractions import Fraction as F
from math import prod,lcm
e=list(combinations(range(4),2))
s=[tuple(int(((m>>i)&1)!=((m>>j)&1)) for i,j in e)for m in range(1,8)]
sites=[[int(i==j)for i in range(6)]for j in range(6)]+s
def basis(A):
 M=[list(map(F,a))+[F(i==j)for j in range(6)]for i,a in enumerate(A)]
 for k in range(6):
  u=next((i for i in range(k,6)if M[i][k]),None)
  if u is None:return
  M[k],M[u]=M[u],M[k]
  v=M[k][k];M[k]=[x/v for x in M[k]]
  for j in range(6):
   if j==k:continue
   q=M[j][k]
   if q:M[j]=[a-q*b for a,b in zip(M[j],M[k])]
 return [a[6:]for a in M]
B=[(I,H)for I in combinations(range(13),6) if (H:=basis([sites[i]for i in I]))]
TT=[]
for I in combinations(range(6),3):
 c={0}
 while 1:
  D=c|{b if a in c else a for i in I for a,b in [e[i]]if a in c or b in c}
  if c==D:break
  c=D
 if len(c)==4:TT.append(I)
assert len(TT)==16
den=lcm(*(q.denominator for _,m in B for row in m for q in row))
assert den<=12
groups={}
for j in range(1,8):
 desc=np.array([[sum(int(den*q)for q,i in zip(row,inds)if i>=6 and not ((i-5)&j).bit_count()%2)
                 for row in h]for inds,h in B],dtype='int64')
 prim=np.array([[sum(int(den*q)for q,i in zip(row,inds)if i>=6 and ((i-5)&j).bit_count()%2)
                 for row in h]for inds,h in B],dtype='int64')
 groups[j]=desc,prim,np.array([(den if ((t&j).bit_count()&1) else 0) for t in range(1,8)],dtype='int64')
cuts=np.array(s,dtype='int64').T
tt=np.array(TT,dtype='int64')
D=100000
scale=den*D
assert all(np.max(abs(x))<=12 for pair in groups.values() for x in pair[:2])
def det(l,h):
 best=None
 for j,(A,B,M) in groups.items():
  pts=A*l+B*h
  ww=pts[(pts>=0).all(axis=1)]
  vv=ww[((ww@cuts)>=(den*l+(h-l)*M)[None,:]).all(axis=1)]
  assert vv.size
  assert vv.max()<5000000
  weights=np.prod(vv[:,tt].astype(object),axis=2).sum(axis=1)
  val=int(min(weights))
  if best is None or val<best[0]:best=val,j,vv[int(np.argmin(weights))]
 return F(best[0],scale**3),best[1:]
def verify(loss, Lbound, witness):
    start=F(1,2)-loss
    end=F(3,2)-loss
    a0,b0=int(start*1000),int(end*1000)
    best=None
    for X in range(a0,b0):
        a,b=F(X,1000),F(X+1,1000)
        A=max(2-b-loss,3*a/4)
        first=F(5,2)-b if b<=F(1,2) else F(9,4)-b/4-1/(16*b)
        H=min(first,3-b)-loss
        assert 0<A<=H
        floor=lambda u:D*u.numerator//u.denominator
        K,att=det(floor(A),floor(H))
        value=a*K
        if best is None or value<best[0]:best=value,X,att[0]
    assert best==witness and best[0]>Lbound
    print('PASS loss',loss, 'min certified vertical Gram determinant',best)
verify(F(0),F(98,100),(F(988455983303,10**12),1142,1))
verify(F(11,100),F(80,100),(F(1635322689,2000000000),1079,1))
# Pinned 4D occupied bound <0.63, pi>3.14.
c=F(2)*F(63,100)/F(314,100)**2
target=F(3**32,5*2**40)
for loss,cap,lower in [(F(0),F(1636,100),F(98,100)),(F(11,100),F(1832,100),F(80,100))]:
    assert 512*c/(2-loss)**2<cap
    assert F(256,1)/lower<cap**2
    assert F(32,1)/(F(3,2)-loss)**2<cap
    assert cap**2<target
assert 256*c/F(135,100)**2<18 and 18**2<target
assert 11*F(5,1000)+4*F(5,1000)**2+12*F(4,1000)<=F(11,100)
assert F(11,4)**2>F(15,2)
print('PASS exact 39D scoped center density caps 16.36 and 18.32')
