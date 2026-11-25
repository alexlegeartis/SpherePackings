"""Exact checker, 34D two-coset A6 + H28 packing. Python 3 and numpy."""
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
B=D["generator_syndromes"];u=D["offset_syndrome"];assert D["half_shift_blocks"]==[0,1,1,1] and D["half_radical_mask"]==82
def tag(m):
 s=0
 for i in range(7):
  if (m>>i)&1:s^=B[i]
 return s
z=costs(0,0);assert (z==0).sum()==1 and min(z[z>0])>=32
assert all(M[i^j]==M[i]^M[j] for i in range(128) for j in range(128))
assert len(B)==7 and D["local_shift"]=={"numerator":3,"denominator":7}
counts={}
for m in range(128):
 w=m.bit_count(); odd=w%2; 
 n=int(min(costs(tag(m)^(u if odd else 0),odd))) # 4 times external norm squared
 loc7=7*(w+2*max(0,3-w))-9 if odd else 7*w
 if m: assert 7*n+4*loc7>=224,(m,n,loc7)
 counts[(w,n,loc7)]=counts.get((w,n,loc7),0)+1
print("nonzero H kernel min:",min(32,min(z[z>0])),"/4")
print("weight count exterior/4 local/7")
for (w,n,local),cnt in sorted(counts.items()):
 print(w,cnt,n,local)
print("minimum=8. covol=32768*sqrt(7); two cosets; center density=8/sqrt(7)")
