#!/usr/bin/env python3
"""Replay the frozen polynomial sources, height matrix, basis and exact cusp glue.
Requires python-flint 0.9.0. No enumeration, reduction or random choices.
"""
from pathlib import Path
from flint import fmpz_mat as Z, fmpq_mat as Q, fmpq as Fr, fq_default_ctx, fq_default_poly_ctx, fmpz_mod_poly_ctx
import gzip,json,hashlib
D=Path(__file__).with_name("gram96.json.gz")
with gzip.open(D,'rt') as s:D0=json.load(s)
p,n,k=47,92,4
K=fq_default_ctx(modulus=fmpz_mod_poly_ctx(p)([1,1,1]),fq_type='FQ_NMOD')
R=fq_default_poly_ctx(K)
def el(i):
 assert isinstance(i,int) and 0<=i<p*p
 return K([i%p,i//p])
t=K([3,23]);w=K([0,1]);c=t**24
assert c==K([4,8]) and t**48==-K(1)
fp={str(c*(a+w*b)):(a,b) for a in range(p) for b in range(p)}
E=[]
labels=[]
for (A,B,C,E0),(xr,yr) in zip(D0['source_parameters'],D0['source_point_polynomials'],strict=True):
 x,y=R(list(map(el,xr))),R(list(map(el,yr)))
 assert x==-R([el(B),el(A)])**16 and y==R([el(E0),el(C)])**24
 assert (not x.is_zero()) and (not y.is_zero())
 assert x.degree()<=16 and y.degree()<=24
 assert y*y-x**3==R([1]+[0]*47+[1])
 E.append((x,y))
 labels.append(fp[str(x(t)/y(t))])
assert len(E)==len(set((str(x),str(y)) for x,y in E))==210
m=len(E);H=Z(m,m)
for i,(x,y) in enumerate(E):
 H[i,i]=16
 for j,(u,v) in enumerate(E[:i]):
  xx,yy=x-u,y-v
  if xx.is_zero() and y==-v:
   hij=-16   # opposite section: height pairing (P,-P)=-16
  else:
   assert not (xx.is_zero() and yy.is_zero())
   f=xx.gcd(yy).degree()
   at_inf=[]
   if not xx.is_zero():at_inf.append(16-xx.degree())
   if not yy.is_zero():at_inf.append(24-yy.degree())
   assert f>=0 and min(at_inf)>=0
   hij=8-f-min(at_inf)
  H[i,j]=H[j,i]=hij
L=Z(D0['horizontal_gram']); V=Z(D0['source_coordinates_in_horizontal_basis'])
B=Q([[Fr(e)for e in row]for row in D0['horizontal_basis_from_independent']])
I=D0['source_indices_independent']
assert len(I)==len(set(I))==n and all(0<=i<m for i in I)
F=Q([[H[i,j]for j in I]for i in I])
assert F.det()>0
assert B*F*B.transpose()==Q(L)
assert Q(V)*B*F==Q([[H[i,j]for j in I]for i in range(m)])
assert V*L*V.transpose()==H and L.det()==p**14
# Compare image of the integer span to the identity row HNF.
h=V.hnf()
assert all(h[i,j]==(1 if i==j else 0) for i in range(m)for j in range(n))
A=Z(D0['cusp_map']);psi=V*A
assert all(psi[i,j]%p==labels[i][j]for i in range(m)for j in range(2))
T=Z(D0['glue_basis_from_horizontal_and_D4'])
D4=Z([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]])
assert D4.det()==4
N=Z(n+k,n+k)
for i in range(n):
 for j in range(n):N[i,j]=L[i,j]
for i in range(k):
 for j in range(k):N[n+i,n+j]=D4[i,j]
cols=((1,0),(26,33),(36,35),(35,37))
C=Z(A.tolist()+[[-a,-b]for a,b in cols])
assert T.det()==p**2 and all(x%p==0 for x in (T*C).entries())
M=Z(D0['gram'])
assert M==T*N*T.transpose() and M.det()==4*p**18
# Digest of canonical full integer 96D Gram, independent of the gzip timestamp.
canonical=json.dumps(D0['gram'],separators=(',',':')).encode()
print('PASS',m,'sections, rank',n+k,'Gram determinant',M.det())
print('Canonical Gram SHA256',hashlib.sha256(canonical).hexdigest())
