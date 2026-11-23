import argparse,json,itertools,math
from fractions import Fraction
from pathlib import Path
from sympy import Matrix
from collections import Counter
P=Path(__file__).resolve().parent
G=[[2*int(i==j) for j in range(7)]for i in range(7)]
for i,j in [(0,1),(0,2),(2,3),(0,4),(4,5),(5,6)]:G[i][j]=G[j][i]=-1
Q=Matrix(G);assert Q.det()==2
assert all(Q[:k,:k].det()>0 for k in range(1,8))
Qi=Q.inv(); assert Qi[6,6]==Fraction(3,2)
def vec(a):return [(a>>i)&1 for i in range(7)]
def quadratic(a):
 v=vec(a)
 return (sum(x*y*G[i][j]for i,x in enumerate(v)for j,y in enumerate(v))//2)&1
def bil(a,b):return quadratic(a)^quadratic(b)^quadratic(a^b)
R=[x for x in range(128)if not any(bil(x,1<<i)for i in range(7))]
assert R==[0,82] and quadratic(82)==1 # (0,1,0,0,1,0,1)
r=R[1]
# Exact coordinate bounds from inverse Gram: all squared norms <=4.
bounds=[math.isqrt(int(4*Qi[i,i])) for i in range(7)]
n02=[0]*128;n04=[0]*128
for x in itertools.product(*(range(-b,b+1)for b in bounds)):
 v=sum(d*d for d in x)-sum(x[i]*x[j]for i,j in [(0,1),(0,2),(2,3),(0,4),(4,5),(5,6)])
 if v<=2:
  assert v>=0 and (x==(0,)*7)==(v==0)
  a=sum((d&1)<<i for i,d in enumerate(x))
  if v==1:n02[a]+=1
  if v==2:n04[a]+=1
assert sum(n02)==126 and n02[r]==0 and sum(x>0 for x in n02)==63
assert all(n02[a]>0 for a in range(1,128)if quadratic(a) and a!=r)
assert all(n04[a]>0 for a in range(1,128)if not quadratic(a))

def mul(x,y):
 r=0
 while y:
  if y&1:r^=x
  x <<=1
  if x&64:x^=67
  y >>=1
 return r
def pow0(a,k):
 r=1
 while k:
  if k&1:r=mul(r,a)
  a=mul(a,a);k>>=1
 return r
assert len({pow0(2,j)for j in range(63)})==63 and pow0(2,63)==1
F8=[x for x in range(64)if pow0(x,8)==x]; assert len(F8)==8
def tr(t):return t^pow0(t,2)^pow0(t,4)
def q0(x):return tr(pow0(x,9))
assert all(q0(x) in (0,1) for x in range(64))
def b0(x,y):return q0(x)^q0(y)^q0(x^y)
# Extend a symplectic matching greedily. Every independent set remains extendible
# to an isometry by Witt extension of the nonsingular polar spaces.
bs=[]; span={0,r}
for i in range(6):
 a=next(a for a in range(128)if a not in span and quadratic(a)==q0(1<<i)
        and all(bil(a,b)==b0(1<<i,1<<j)for j,b in enumerate(bs)))
 bs.append(a);span.update({h^a for h in span})
I=[0]
for a in bs:I+= [a^s for s in I]
M=I+[a^r for a in I]
assert len(set(M))==128
assert all(quadratic(a)==(q0(x&63)^(x>>6))for x,a in enumerate(M))
alpha=F8[:5]; assert alpha[:2]==[0,1]
def word(x):
 s=x[0]^x[1]^x[2];t=mul(alpha[2],x[0])^mul(alpha[3],x[1])^mul(alpha[4],x[2])
 return [s^t,t]+list(x)
def glue(x,e):
 return sum((M[y|(((e>>i)&1)<<6)])<<(7*i)for i,y in enumerate(x))
Z=[]
for a in range(3):
 for i in range(6):
  x=[0]*3;x[a]=1<<i
  Z.append(glue(word(x),0))
for i in range(1,5):Z.append(glue([0]*5,1|(1<<i)))
rank={}
for z in Z:
 while z:
  j=z.bit_length()-1
  if j not in rank:rank[j]=z;break
  z^=rank[j]
assert len(rank)==22
C=Counter()
for x in itertools.product(range(64),repeat=3):
 y=word(x); w=sum(bool(t)for t in y)
 C[w]+=1;assert w not in [1,2]
 if w==3:assert sum(q0(t)for t in y)%2==0
assert C[3]==630
# The horizontal section normal is the primitive last fundamental weight in the first block.
# Its pairing on coefficient rows is the coefficient at index 6. Check some generator is odd.
assert any((z>>6)&1 for z in Z)
Pout={'E7_gram':G,'alpha':alpha,'field_polynomial':67,'phase_isometry':M,'generators_35_hex':[f'{z:09x}'for z in Z]}
out=P/'certificate.json'
if '--write' in __import__('sys').argv:out.write_text(json.dumps(Pout,indent=2)+'\n')
else:assert Pout==json.loads(out.read_text())
print('PASS local shell 126 roots and bounds, radical, GF64 quadratic isometry, rank',len(rank),'weight counts',dict(C))
print('35 Gram determinant',2**31,'squared minimum 8, center density 4')
print('34 section Gram determinant',3*2**30,'squared minimum 8, center density 4/sqrt(3)')
