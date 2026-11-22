from pathlib import Path
from itertools import combinations
import json,re
import numpy as np
from sympy import Matrix
HERE=Path(__file__).resolve().parent
c=json.loads((HERE/'certificate.json').read_text())
manifest={a['dimension']:a for a in json.loads((HERE/'numerators.json').read_text())}
results={a['dimension']:a for a in json.loads((HERE/'results.json').read_text())['rows']}
import hashlib
S=HERE.parent/'novel_routes_20261113'/'higher_parent'/'neb72.html'
assert hashlib.sha256(S.read_bytes()).hexdigest()=='2734b670349529e98de95af8810867077eb8d802363731b59cdb8082c582b9fa'
J=np.array(c['leech_frame'],dtype='i8')
rawtext=S.read_text().split('NAME="GRAM"',1)[1].split('NAME="DET"',1)[0]
n,*vals=map(int,re.findall(r'-?\d+',re.sub('<[^>]*>',' ',rawtext)))
assert n==72 and len(vals)==n*n
GG=Matrix(n,n,vals)
assert GG.det(method='domain-ge')==1 and GG.T==GG and all(GG[i,i]%2==0 for i in range(72))
G=np.array(GG.tolist(),dtype='i8')
def exact(x):
    assert all(z.denominator==1 for row in x.tolist() for z in row)
    return Matrix([[int(z) for z in row] for row in x.tolist()])
U=Matrix(c['leech_basis_change'])
Ui=exact(U.inv(method="DM")); assert abs(U.det())==1
I=Matrix([[int(i==j) for j in range(24)] for i in range(24)])
F0=exact(Matrix(G[:24,:24].tolist())/2)
A0=exact(Matrix(G[:24,48:].tolist())*F0.inv(method="DM"))
B0=I-A0
F=np.array((U*F0*U.transpose()).tolist(),dtype='i8')
A=np.array((U*A0*Ui).tolist(),dtype='i8')
B=np.eye(24,dtype='i8')-A
O=np.zeros((24,24),'i8')
ii=np.array(Ui.tolist(),dtype='i8')
t=np.array((A0*Ui).tolist(),dtype='i8')
b=np.array((B0*Ui).tolist(),dtype='i8')
T=np.block([[ii,ii,t],[O,b,b],[O,O,2*ii]])
Ti=Matrix(T.tolist()).inv(method="DM")
assert all(a.denominator<=2 for row in Ti.tolist() for a in row)
# The Leech block F is as in Nebe's exact construction (first G block divided by 2).
# Index of T is 2^36 because T diag(F,F,F) T^T = 2 G and both G and F are unimodular.
blockF=np.kron(np.eye(3,dtype='i8'),F)
assert np.array_equal(T@blockF@T.T,2*G)
assert Matrix(F.tolist()).det()==1
assert np.array_equal(J@F@J.T,8*np.eye(24,dtype='i8'))
A%=2; B%=2
I=np.eye(24,dtype='i8')
assert np.array_equal((A+B)%2,I) and np.array_equal(A@A%2,A)
assert (A@F@A.T%2==0).all() and (B@F@B.T%2==0).all()
assert ((A@F@A.T).diagonal()%4==0).all() and ((B@F@B.T).diagonal()%4==0).all()
p24=1<<np.arange(24)
def wpack(v):return int(v%2@p24)
def rank(words):
 local={}
 for x in words:
  while x:
   i=x.bit_length()-1
   if i not in local:local[i]=x;break
   x^=local[i]
 return len(local)
assert rank(map(wpack,A))==rank(map(wpack,B))==12
for x,y,z in T.reshape(72,3,24):
 assert np.array_equal(x@A%2,y@A%2) and np.array_equal(x@A%2,z@A%2)
 assert ((x+y+z)@B%2==0).all()
# These three independent rank-12 constraints characterize the entire T image.
# Enumerate the complete norm-four Leech shell (integer vectors in a frame, sum of squares 32).
fwd=F@J.T
p2=np.unique(fwd%2,axis=0)
assert [r.tolist() for r in p2]==[[0]*24,[1]*24]
def pack4(r): return (wpack(r),wpack(r>>1))
def add4(a,b): return (a[0]^b[0],a[1]^b[1]^(a[0]&b[0]))
R={(0,0)}
for r in fwd:
 a=pack4(r%4); a2=add4(a,a)
 R={add4(t,x) for t in R for x in [(0,0),a,a2,add4(a2,a)]}
assert len(R)==8192 and all(lo or hi.bit_count()!=4 for lo,hi in R)
rows=[]
for i,j in combinations(range(24),2):
 for x,y in [(4,4),(4,-4),(-4,4),(-4,-4)]:
  q=np.zeros(24,dtype='i8');q[i]=x;q[j]=y;rows.append(q)
for lo,hi in R:
 if lo==0 and hi.bit_count()==8:
  idx=[i for i in range(24) if hi>>i&1]
  q=np.zeros(24,dtype='i8')
  for z in range(256):
   q[idx]=[-2 if z>>k&1 else 2 for k in range(8)]
   rows.append(q.copy())
 elif lo==(1<<24)-1:
  q=np.array([-1 if hi>>i&1 else 1 for i in range(24)],dtype='i8')
  for j in range(24):
   v=q.copy();v[j]*=-3;rows.append(v)
v=np.array(rows)@J
W=v[(v%8==0).all(axis=1)]//8
assert len(W)==len(np.unique(W,axis=0))==196560 and ((((W@F)*W).sum(1))==4).all()
M4={}
for a,b in zip((W@A%2)@p24,(W@B%2)@p24):
 M4.setdefault(int(a),set()).add(int(b))
assert len(M4)==4095 and 0 not in M4 and set(map(len,M4.values()))=={24}
def bits(i):return (i>>np.arange(24))&1
def orth(a,b):return int(bits(a)@F@bits(b))%2==0
def syndrome(r):
 x,y,z=np.array(r,dtype='i8').reshape(3,24)
 return wpack((x+z)@A),wpack((y+z)@A),wpack((x+y+z)@B)
def prove(r,cut):
 a,b,t=syndrome(r)
 if cut==4:
  assert a==b==0 and t!=0
  # A=0 has no norm-four label, and B is singular, so no single norm-six either.
  return
 assert cut==6 and a and b and a!=b
 U0=(a,b,a^b)
 assert all(orth(s,t) for s in U0)
 assert all(orth(s,x) for j in U0 for x in M4[j] for s in U0)
 # For two occupied blocks their two A labels are any of (a,b),(a,a+b),(b,a+b).
 assert all(t!=x^y for a0,b0 in combinations(U0,2) for x in M4[a0] for y in M4[b0])
 # This excludes 4+4. For 4+6, the norm-four B label and the total B label
 # annihilate U0, so the norm-six B label annihilates its own A label, impossible.
def build(raw,vertical):
 r=len(raw)
 H=np.array(raw,dtype='i8')
 assert (H==H%2).all()
 w=2*Matrix(H.tolist())*Ti
 W=exact(w)
 vg=Matrix(vertical)
 X=W*GG/2
 N=W*GG*W.transpose()/4+vg
 qq=Matrix([v+x for v,x in zip(GG.tolist(),X.transpose().tolist())]+[v+x for v,x in zip(X.tolist(),N.tolist())])
 assert qq.det(method='domain-ge')==vg.det()
 det=vg.det()
 row=manifest[72+r]
 assert row['vertical_gram']==vertical and row['numerator']==W.tolist() and row['determinant']==str(det)
 print('PASS',72+r,'full lattice, min squared 8, Gram determinant',det,flush=True)
prove(c['raw73'],6)
build([c['raw73']],[[2]])

from itertools import product
from math import isqrt
vertical=[[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]]
V0=Matrix(vertical)
rootpar={tuple(x%2 for x in z) for z in product(range(-3,4),repeat=4)
         if (Matrix([z])*V0*Matrix(z))[0]==2}
nonroot=set(product(range(2),repeat=4))-rootpar-{(0,0,0,0)}
assert len(rootpar)==12 and len(nonroot)==3
ks=sorted(nonroot)
basis=next([a,b,*ks[:2]] for a,b in combinations(sorted(rootpar),2)
           if Matrix([a,b,*ks[:2]]).det()%2)
q=Matrix(basis).inv_mod(2)
rawD4=np.array(q.tolist(),dtype='i8')@np.array(c['raw_D4_seed'],dtype='i8')%2

def verify_short(raw,vertical):
 r=len(raw)
 V=Matrix(vertical)
 Vi=V.inv(); assert V.is_positive_definite
 bounds=[isqrt((8*Vi[i,i]).p//(8*Vi[i,i]).q) for i in range(r)]
 checked=set()
 for z in product(*(range(-b,b+1) for b in bounds)):
  if not any(z):continue
  zz=np.array(z,dtype='i8')
  v=int(zz@np.array(vertical,dtype='i8')@zz)
  if v>=8:continue
  assert v in (2,4,6)
  raw0=zz@np.array(raw,dtype='i8')%2
  if (tuple(raw0),v) in checked:continue
  checked.add((tuple(raw0),v))
  a,b,t=syndrome(raw0)
  if v==2:
   prove(raw0,6)
  elif v==4:
   if a==b==0:
    prove(raw0,4)
   elif a==b or not(a and b):
    assert t not in M4[a or b] and orth(a or b,t)
   # (a,b) independent prohibits any single-block representation.
  else:
   assert (a,b,t)!=(0,0,0) # one occupied block has squared raw norm ≥4
 print('vertical short parities',r,len(checked))

for r in (2,3,4):
 raw=rawD4[:r]
 vertical_r=np.array(vertical)[:r,:r].tolist()
 verify_short(raw,vertical_r)
 build(raw.tolist(),vertical_r)

raw0=np.array(c['raw73'],dtype='i8').reshape(3,24)
Y=raw0[[1,0,2]]
A2=np.vstack([Y.ravel(),np.roll(Y,1,axis=0).ravel()])
full=np.vstack([rawD4,A2])
V6=np.zeros((6,6),'i8')
V6[:4,:4]=vertical
V6[4:,4:]=[[2,-1],[-1,2]]
for r in (5,6):
  vv=V6[:r,:r].tolist()
  verify_short(full[:r],vv)
  build(full[:r].tolist(),vv)

from fractions import Fraction as Q
from math import factorial, prod, isqrt
def atan(q,n):
  return sum(Q((-1)**i,(2*i+1)*q**(2*i+1)) for i in range(n+1))
pi=16*atan(5,41)-4*atan(239,40)
assert pi>3
for n,a in results.items():
  v=2**n/Q(a['gram_determinant'])
  assert manifest[n]['determinant']==a['gram_determinant']
  assert v==Q(a['center_density_squared'])
  ball=(pi**n/factorial(n//2)**2 if n%2==0 else 2**(n+1)*pi**(n-1)/prod(range(1,n+1,2))**2)
  assert Q(a['occupied_density_floor'])**2<ball*v
  for which in ('unrestricted','lattice'):
    p=Q(a[f'prior_{which}_center_density_squared'])
    pct=Q(a[f'gain_over_prior_{which}_percent_floor'])
    assert (1+pct/100)**2<=v/p<(1+(pct+Q(1,10**6))/100)**2
  ext=Q(a['named_external_comparison_center_density_squared'])
  pct=Q(a['gain_over_named_external_comparison_percent_floor'])
  assert (1+pct/100)**2<=v/ext<(1+(pct+Q(1,10**6))/100)**2
  print('PASS exact density comparisons',n)
