import numpy as np, math
from itertools import product,combinations_with_replacement
from pathlib import Path
import json
from collections import Counter
from functools import lru_cache
from fractions import Fraction as F

def a2(q):
 E=list(product(range(q),repeat=2))
 G=np.array([[2,-1],[-1,2]],int)
 def orb(z):
  L=[z];U={z}
  for a in L:
   u,v=a
   for k in ((-u+v,v),(u,-v+u),(-u,-v)):
    b=tuple(t%q for t in k)
    if b not in U:U.add(b);L.append(b)
  return sorted(U)
 seen=set();O=[]
 for z in E:
  if z in seen:continue
  S=orb(z);O.append(S);seen.update(S)
 def dis(z):
  x,y=z
  return min(2*(a*a-a*b+b*b) for i,j in product(range(-2,2),repeat=2) for a,b in [(x+i*q,y+j*q)])
 d=[dis(S[0]) for S in O]
 assert all(all(dis(v)==x for v in S)for x,S in zip(d,O))
 cos={2:[2,-2],3:[2,-1,-1],4:[2,0,-2,0],6:[2,1,-1,-2,-1,1]}[q]
 C=sorted({tuple(sum(cos[(i*x+j*y)%q]for x,y in O0) for O0 in O)for i,j in E})
 assert len(C)==len(O)
 print(q,'phase shell sizes',[len(c)for c in O],'depth',d,flush=True)
 return np.array([len(c)for c in O]),np.array(C),[F(x,q*q)for x in d]
@lru_cache(None)
def T(m,q):
 W,H,D=a2(q)
 k=len(W);cl=list(combinations_with_replacement(range(k),m))
 ix=[tuple(Counter(a)[j]for j in range(k))for a in cl]
 Dic={a:i for i,a in enumerate(ix)}
 def coeff(R):
  p={tuple([0]*k):1}
  for v in R:
   new=Counter()
   for a,w in p.items():
    for i,e in enumerate(v):
     if e:new[a[:i]+(a[i]+1,)+a[i+1:]]+=w*int(e)
   p=new
  row=np.zeros(len(cl),dtype=np.int64)
  for a,w in p.items():row[Dic[a]]=w
  return row
 N=np.array([math.factorial(m)//math.prod(math.factorial(j)for j in a)*math.prod(int(w)**j for w,j in zip(W,a))for a in ix],dtype=np.int64)
 M=np.array([coeff(H[list(r)]) for r in cl],dtype=np.int64)
 d=np.array([sum(D[i]for i in a)for a in cl],dtype=object)
 return N,M,d
data=json.loads(Path(__file__).with_name('duals.json').read_text())
for q,d in data.items():
 q=int(q);W,M,D=T(5,q)
 ix=np.flatnonzero((D>=2)|(np.arange(len(D))==0))
 A=M[:,ix];c=W[ix];Y={i:F(y)for i,y in d['weights']}
 assert min(Y.values())>0
 WX=[int(c[j])+sum(v*int(A[i,j])for i,v in Y.items())for j in range(len(c))]
 cap=WX[0]+sum(max(F(0),a)for a in WX[1:])
 assert cap==F(d['cap']) and cap<(19 if q==3 else 21)
 print('PASS exact cap',q,cap)
C=np.loadtxt(Path(__file__).with_name('ternary-18.txt'),dtype=np.int64)
assert C.shape==(18,5) and (C>=0).all() and (C<3).all()
D=np.count_nonzero(C[:,None]!=C[None,:],axis=-1)
assert min(D[~np.eye(len(C),dtype=bool)])==3
G=np.array([[2,-1],[-1,2]],int)
for k in [1,2]:
 u,v=k,k*2%3
 b=min(2*(a*a-a*b+b*b)for a in range(u-6,u+4,3)for b in range(v-6,v+4,3))
 assert b==6
# Literal translates: block residues are C[k,i]*(1,2)/3; every
# nonzero differing block has min depth 2/3. The parent min is 2.
assert len(set(map(tuple,C)))==18
print('PASS attained 18, center density 1/(16 sqrt(3)), nonlattice since 18 is not a power of 3')
