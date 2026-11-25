"""Exact checker, 34D D6/E7 index-two lattice. Python 3 and numpy."""
import itertools, json, math
from pathlib import Path
import numpy as np
D=json.loads(Path(__file__).with_name("cert.json").read_text())
G=np.array(D["E7_gram"],dtype=int); M=D["phase_isometry"]; A=D["alpha"][1:]
F=np.array([[2,-1,-1,0,-1,0,0],[0,3,-1,0,-1,0,0],[0,0,4,-3,-2,0,0],[0,0,0,5,-2,0,0],[0,0,0,0,4,-5,0],[0,0,0,0,0,3,-4],[0,0,0,0,0,0,1]])
c=[30,10,5,3,3,5,40];assert np.array_equal(F.T@np.diag(c)@F,60*G)
def det(a):
 if len(a)==1:return int(a[0][0])
 return sum((-1)**i*int(t)*det(np.delete(a[1:],i,axis=1)) for i,t in enumerate(a[0]) if t)
assert det(G)==2 and sorted(M)==list(range(128)) and M[64]==82 and D["field_polynomial"]==67
def gf(x,y):
 r=0
 while y:
  if y&1:r^=x
  x*=2
  if x&64:x^=67
  y//=2
 return r
t=np.array([[gf(i,j) for j in range(64)] for i in range(64)])
C=[{},{}];x=[0]*7;r=[(82>>i)&1 for i in range(7)]
def enum(j,budget):
 if j<0:
  n=(1440-budget)//60
  for m in (0,1):
   if any((v-m*z)%2 for v,z in zip(x,r)):continue
   k=sum((((v-m*z)//2)%2)<<i for i,(v,z) in enumerate(zip(x,r)))
   C[m][k]=min(n,C[m].get(k,999))
  return
 off=sum(int(F[j,i])*x[i] for i in range(j+1,7));d=int(F[j,j]);h=math.isqrt(budget//c[j])
 for y in range(-((h+off)//d),(h-off)//d+1):
  x[j]=y;enum(j-1,budget-c[j]*(d*y+off)**2)
enum(6,1440)
assert [len(v) for v in C]==[128,128]; L=np.array([[v[k] for k in M] for v in C]); assert max(L[0])<=24 and max(L[1])<=24
j=np.arange(32768); k0=j&63;k1=(j>>6)&63
h0=(j>>12)&1;h1=(j>>13)&1;h2=(j>>14)&1
inv=next(j for j in range(64) if gf(A[2]^A[3],j)==1)
def costs(s,typ):
 a=(s&63)^k0^k1; b=((s>>6)&63)^t[A[0],k0]^t[A[1],k1];k2=t[inv,b^t[A[3],a]];k3=a^k2
 h3=(s>>12)^h0^h1^h2
 return L[0,k0+64*h0]+L[typ,k1+64*h1]+L[typ,k2+64*h2]+L[typ,k3+64*h3]
B=D["phase_syndromes"];u=D["offset_syndrome"];assert D["local_lattice"]=="D6" and len(B)==6 and D["half_shift_blocks"]==[0,1,1,1] and D["half_radical_mask"]==82
def tag(x):
 return __import__('functools').reduce(int.__xor__,(B[i] for i in range(6) if x>>i&1),0)
z=costs(0,0);assert (z==0).sum()==1 and min(z[z>0])>=32
assert all(M[i^j]==M[i]^M[j] for i in range(128) for j in range(128))
assert tag(63)==4096
# Enumerate integer local coordinates near the origin independently in original D6 coordinates.
# This box contains every norm-minimizer needed: squared local distances at most six
# cannot have an integral coordinate of magnitude >2 or half coordinate of magnitude >5/2.
local=[[100]*64 for j in (0,1)]
for y in itertools.product(range(-2,4),repeat=6):
 s=sum(y)
 if s%2: continue
 x=sum((v%2)<<i for i,v in enumerate(y[:5]))|(((s//2)%2)<<5)
 for m in (0,1):
  d=sum((2*v-m)**2 for v in y)
  local[m][x]=min(local[m][x],d)
hist={}
for x in range(64):
 for m in (0,1):
  loc=local[m][x];assert loc<=24
  n=int(min(costs(tag(x)^(u if m else 0),m)))
  if x or m:assert n+loc>=32,(m,x,n,loc)
  key=(m,loc,n);hist[key]=hist.get(key,0)+1
print("nonzero H kernel min",min(32,min(z[z>0])),"/4; 2t phase",tag(63))
print("bank type, local norm times 4, exterior norm times 4, count")
for key,v in sorted(hist.items()):print(*key,v)
print("minimum = 8; full lattice covolume = 32768; center density = 4")
