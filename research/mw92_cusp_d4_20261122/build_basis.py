from pathlib import Path
from collections import defaultdict
import json,sys,time,numpy as np
from flint import fmpz_mat as Z, fmpq_mat as Q, nmod_mat, fq_default_ctx, fq_default_poly_ctx, fmpz_mod_poly_ctx
p=47; rank=92; den=14; Chi=8; q0=p*p
field=fq_default_ctx(modulus=fmpz_mod_poly_ctx(p)([1,1,1]),fq_type='FQ_NMOD')
def x(i):return field([int(i)%p,int(i)//p])
def enc(f):
 s=f.to_list()
 return (int(s[0]) if s else 0)+(int(s[1])*p if len(s)>1 else 0)
ring=fq_default_poly_ctx(field)
A=[x(i) for i in range(q0)]
G=defaultdict(list)
for i,a in enumerate(A):G[enc((a**(p+1)))].append(i)
rng=np.random.default_rng(38)
nb=int(sys.argv[1])if len(sys.argv)>1 else 210
rows=[];seen=set();poly=[]
for _ in range(nb):
  while True:
   a=int(rng.integers(q0)); n0=enc(A[a]**(p+1));b=int(rng.choice(G[(1-n0)%p]))
   e=int(rng.choice(G[1]));c=enc((-(A[b]**p))*A[e]);d=enc((A[a]**p)*A[e])
   aa,bb,cc,dd=map(x,[a,b,c,d])
   assert aa**(p+1)+bb**(p+1)==1 and cc**(p+1)+dd**(p+1)==1 and aa*cc**p+bb*dd**p==0
   xx=-(ring([bb,aa])**16); yy=ring([dd,cc])**24
   # Exact duplicate polynomials are not needed.
   tag=(str(xx),str(yy))
   if tag not in seen:
    seen.add(tag)
    assert yy*yy==xx**3+ring([1]+[0]*47+[1])
    rows.append([a,b,c,d]);poly.append((xx,yy));break
H=np.zeros((nb,nb),dtype='int64')
for i,(x0,y0) in enumerate(poly):
 H[i,i]=16
 for j,(x1,y1) in enumerate(poly[:i]):
  xx,yy=x0-x1,y0-y1
  if xx.is_zero() and (y0+y1).is_zero():v=-16
  else:
   inter=xx.gcd(yy).degree()
   hi=[]
   if not xx.is_zero():hi.append(16-xx.degree())
   if not yy.is_zero():hi.append(24-yy.degree())
   assert inter>=0 and hi and min(hi)>=0
   v=Chi-inter-min(hi)
  H[i,j]=H[j,i]=v
print('computed height matrix',nb,flush=True)
a=nmod_mat(H.tolist(),101)
# Incremental basis of rows modulo 101.
D0={};good=[]
for j,row in enumerate(H):
 u=row.copy()%101
 while np.any(u):
  k=np.flatnonzero(u)[0]
  if k in D0:u=(u-u[k]*D0[k])%101
  else:D0[k]=(u*pow(int(u[k]),99,101))%101;good.append(j);break
 if len(good)==rank:break
assert len(good)==rank
F0=Z(H[np.ix_(good,good)].tolist());assert F0.det()>0
F=Q(F0)
Y=Q(H[:,good].tolist())*F.inv()
d=1
from math import lcm
for y in Y.entries():d=lcm(d,int(y.denominator))
YI=Z([[int(x*d)for x in r]for r in Y.tolist()])
h=YI.hnf()
P=Q(h.tolist()[:rank])/d;I=P.inv()
P0=P*F*P.transpose();Ys=Y*I
assert all(x.denominator==1 for x in P0.entries()+Ys.entries())
assert P0.det()==p**den
assert Ys*P0*Ys.transpose()==Q(H.tolist())
print('basis source complete; reducing size',flush=True)
P0=Z([[int(x)for x in row]for row in P0.tolist()])
try:
 P0,U=P0.lll(transform=True,delta=.99,rep='gram',gram='approx')
except Exception:
 P0,U=P0.lll(transform=True,delta=.99,rep='gram',gram='exact')
assert abs(U.det())==1
Ys=Ys*Q(U).inv()
assert Ys*P0*Ys.transpose()==Q(H.tolist())
# Recover linear cusp evaluation on the verified integral basis over F47.
R0=[]
t=field([3,23]); c0=t**24; w=field([0,1])
assert enc(c0)==4+8*p
Qbasis={enc(c0*(a+w*b)):(a,b) for a in range(p) for b in range(p)}
for X,Y in poly:
 assert Y(t)!=field(0)
 R0.append(Qbasis[enc(X(t)/Y(t))])
V=np.array([[int(x)for x in r]for r in Ys.tolist()],dtype='i8')
R0=np.array(R0,dtype='i8')
D1={}
for j,row in enumerate(V):
 u=np.r_[row,R0[j]]%p
 while np.any(u[:rank]):
  k=np.flatnonzero(u[:rank])[0]
  if k in D1:u=(u-u[k]*D1[k])%p
  else:D1[k]=(u*pow(int(u[k]),-1,p))%p;break
 if len(D1)==rank:break
assert len(D1)==rank
R=np.zeros((rank,2),dtype='i8')
for j in range(rank-1,-1,-1):
 R[j]=(D1[j][-2:]-D1[j][j+1:rank]@R[j+1:])%p
assert np.array_equal(V@R%p,R0)
print('cusp map verified in working basis; tested original labels',len({tuple(q)for q in R0}),flush=True)
from itertools import combinations
pairs=[(i,j)for i,j in combinations(range(nb),2)if
  (H[i,j]==7 and np.array_equal(R0[i],R0[j])) or
  (H[i,j]==-7 and np.array_equal(R0[i],-R0[j]%p))]
if pairs: print('cusp period witness of norm18 from section pair:',pairs[0],flush=True)
D4=Z([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]])
assert D4.det()==4
Q0=Q(U)
Q=Z(rank+4,rank+4)
for i in range(rank):
 for j in range(rank):Q[i,j]=P0[i,j]
for i in range(4):
 for j in range(4):Q[rank+i,rank+j]=D4[i,j]
C=np.vstack([R,-np.array([[1,0],[26,33],[36,35],[35,37]])])%p
B=Z(rank+4,rank+4)
for i in range(rank+4):
 if i in (rank,rank+1):B[i,i]=p
 else:
  a,b=map(int,C[i]); y=b*pow(33,-1,p)%p
  x=(a-26*y)%p
  B[i,i]=1;B[i,rank]=x;B[i,rank+1]=y
assert B.det()==p**2
assert not np.any(np.array(B.tolist(),dtype='i8')@C%p)
G96=B*Q*B.transpose();assert G96.det()==4*p**18
print('PASS explicit Gram basis rank96, det',G96.det(),flush=True)
import gzip
out={'gram':G96.tolist(),'source_parameters':rows,'horizontal_gram':P0.tolist(),'cusp_map':R.tolist(),
 'source_point_polynomials':[[[enc(e) for e in X.coeffs()],[enc(e) for e in Y.coeffs()]] for X,Y in poly],
 'source_indices_independent':good, 'horizontal_basis_from_independent':[[str(e)for e in row]for row in (Q0*P).tolist()],
 'source_coordinates_in_horizontal_basis':V.tolist(),
 'glue_basis_from_horizontal_and_D4':[[int(e)for e in row]for row in B.tolist()]}

out['gram']=[[int(e)for e in r]for r in out['gram']];out['horizontal_gram']=[[int(e)for e in r]for r in out['horizontal_gram']]
with gzip.open('gram96.json.gz','wt')as f:json.dump(out,f,separators=(',',':'))
