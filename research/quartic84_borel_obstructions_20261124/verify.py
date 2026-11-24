from fractions import Fraction as F
from math import gcd,isqrt
import sys,json
from collections import Counter

from receipt import read
d=read['canonical']
G,R,J,T,V= (d[name]for name in ('G','R','J','U','F'))
def identity(n):return [[int(i==j)for j in range(n)]for i in range(n)]
def neg(A):return [[-x for x in r]for r in A]
def transpose(A):return list(map(list,zip(*A)))
def row(v,A):return [sum(a*b for a,b in zip(v,w))for w in zip(*A)]
def mul(a,b):return [row(v,b)for v in a]
def q(x,G):return sum(a*b for a,b in zip(x,row(x,G)))
def pow(A,t):
 out=identity(len(A))
 while t:
  if t&1:out=mul(out,A)
  t//=2
  if t:A=mul(A,A)
 return out
def det(M,positive=False):
 A=[r.copy()for r in M];n=len(A);D=1;c=1
 for i in range(n-1):
  if not A[i][i]:
   assert not positive
   z=next((k for k in range(i+1,n)if A[k][i]),None)
   if z is None:return 0
   A[i],A[z]=A[z],A[i];c=-c
  p=A[i][i]
  if positive:assert p>0
  for j in range(i+1,n):
   for k in range(i+1,n):
    a=A[j][k]*p-A[j][i]*A[i][k];v,r=divmod(a,D);assert r==0;A[j][k]=v
   A[j][i]=0
  D=p
 if positive:assert A[-1][-1]>0
 return c*A[-1][-1]
Id=identity(84)
assert transpose(G)==G and det(G,True)==2**42
assert pow(R,21)==Id and pow(R,3)!=Id!=pow(R,7)
assert mul(mul(R,G),transpose(R))==G
f=d['known_mod33']['character'];assert all(a==b[0]for a,b in zip(f,mul(R,[[a]for a in f])))
assert gcd(33,*f)==1
for w0 in d['new_half_shell'].values():
 for p,x in w0['trace_zero'].items():
  assert q(x,G)==10
  p=int(p);S=pow(R,21//p);s0=x.copy();trace=x.copy()
  for _ in range(1,p):
   s0=row(s0,S);trace=[a+b for a,b in zip(trace,s0)]
  assert not any(trace)
  print('PASS exact zero order',p,'fixed projection, literal q10')
P=d['known_mod33']['phase_labels'];labels={(a-b)%33 for a in P for b in P}
for g in d['literal_blockers'].values():
 assert set(g['pooled0'])=={f'mod33-{i}'for i in labels}|{'index17-0','index17-1'}
 for k,x in g['pooled0'].items():
  assert q(x,G)==10
  if k.startswith('mod33'):assert sum(a*b for a,b in zip(x,f))%33==int(k[6:])
  else:assert sum(a*b[0]for a,b in zip(x,d['index17_forms'][int(k[-1])]))%17==0
print('PASS literal q10 blockers for all five mod33 differences and both index17 maps')
# The native cyclic antiunitary simplification certifies two half shapes.
assert pow(T,43)==Id!=T and mul(mul(T,G),transpose(T))==G
assert mul(J,J)==neg(Id) and mul(mul(J,G),transpose(J))==G and mul(J,T)==mul(T,J)
assert mul(V,V)==Id and mul(mul(V,G),transpose(V))==G and mul(J,V)==neg(mul(V,J))
assert mul(mul(V,T),V)==pow(T,42)
U=mul(V,J)
for k,A in [('F',V),('FJ',U)]:
 s=d[k+'_half'];C,H=s['C'],s['H'];D=s['D'];N=s['full_primitive_basis']
 assert len(C)==42 and N[-42:]==C and abs(det(N))==1
 assert mul(C,A)==C and mul(D,A)==neg(D) and D==mul(C,J)
 assert mul(mul(C,G),transpose(C))==mul(mul(D,G),transpose(D))==H and det(H,True)==2**63
 assert not any(x for r in mul(mul(C,G),transpose(D))for x in r)
 print('PASS saturated exact 42+42 split for',k)
 if '--enum-input' not in sys.argv:continue
 # Exact floor-LDL bounds. B^{-1}_{ii} from ratios of positive PD minors.
 n=len(H);scl=1<<20;dd=[];L=identity(n);hd=det(H,True)
 for i in range(n):
  c=F(H[i][i])-sum(L[i][j]**2*dd[j]for j in range(i));assert c>0;dd.append(c)
  for j in range(i+1,n):
   L[j][i]=(F(H[j][i])-sum(L[j][a]*L[i][a]*dd[a]for a in range(i)))/c
 W=[isqrt(int(10*F(det([[x for j,x in enumerate(row0)if j!=i]for k,row0 in enumerate(H)if k!=i],True),hd)))for i in range(n)]
 with open(k+'_10.in','w')as o:
  print(n,10,file=o)
  for seq in ([x for r in H for x in r],[int(x*scl)for x in dd],[(x*scl)//1 for r in L for x in r],W):print(*seq,file=o)
for nm,pair0 in d['integral_image'].items():
 x,y=pair0['x'],pair0['y']
 assert q(x,G)==10 and [b-a for a,b in zip(y,row(y,R))]==x
 print('PASS literal integer (R-I) image -> any modulus, source',nm)
for name,a in d['extra_indices'].items():
 p,k,x,f=a['p'],a['a'],a['x'],a['f']
 assert all(v[0]%p==k*b%p for v,b in zip(mul(R,[[b]for b in f]),f))
 assert q(x,G)==10 and sum(b*c for b,c in zip(x,f))%p==0
 print('PASS explicit index',p,'eigencharacter',k,'blocked')
w=read['mutation13'];A=w['source_raw_gram'];x=w['generator_v43'];H=w['H']
assert transpose(A)==A and det(A,True)==43**20*2**42
u=[next(i for i,y in enumerate(r)if y)for r in H]
assert q(x,A)==8*43**2 and all((a-b)%43==0 for a,b in zip(x,row([x[i]%43 for i in u],H)))
print('PASS literal Borel mut. #13 norm8 in its exact Gaussian glue')
