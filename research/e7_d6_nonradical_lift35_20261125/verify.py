"""Exact tables for the 35D lift. Python 3 + numpy. No network. No floating arithmetic."""
import json, math, itertools
from pathlib import Path
from collections import Counter
import numpy as np
d=json.loads(Path(__file__).with_name("cert.json").read_text())
G=np.array(d["E7_gram"],dtype=np.int64)
F=np.array([[2,-1,-1,0,-1,0,0],[0,3,-1,0,-1,0,0],[0,0,4,-3,-2,0,0],[0,0,0,5,-2,0,0],[0,0,0,0,4,-5,0],[0,0,0,0,0,3,-4],[0,0,0,0,0,0,1]],dtype=np.int64)
D=[30,10,5,3,3,5,40]
assert np.array_equal(F.T@np.diag(D)@F,60*G)
M=d["phase_isometry"]; assert sorted(M)==list(range(128)) and M[64]==82 and all(M[i^j]==M[i]^M[j] for i in range(128) for j in range(128))
A=d["alpha"][1:];B=d["phase_syndromes"];of=d["offset_syndrome"]
T=[82*int(i) for i in d["half_shift_blocks"]]; hloc=d["new_local_shift_times_2"]; bIn=d["new_half_masks"];u=d["new_integral_syndrome"]
assert hloc==[3,0,0,0,0,1] and bIn==[21,60,41,40] and u==896 and d['new_layer_height_squared']==[3,2]
def fld(x,y):
 z=0
 while y:
  if y%2:z^=x
  x*=2
  if x&64:x^=d["field_polynomial"]
  y//=2
 return z
mul=np.array([[fld(i,j) for j in range(64)] for i in range(64)])
def syn(j,mask):
 t=M.index(mask);return (t&63)+64*fld(A[j],t&63)+4096*(t>>6)
def fx(x):
 v=0
 for i in range(6):
  if x>>i&1:v^=B[i]
 return v
# The local array is the exact least squared norm of doubled E7 coordinate vectors in each mod-4 residue.
tab=np.full(4**7,1000,dtype=np.int64); V=[0]*7
FF=F.tolist()
def enum(j,budget):
 if j<0:
  idx=sum((v%4)<<(2*i) for i,v in enumerate(V)); tab[idx]=min(int(tab[idx]),(1440-budget)//60)
  return
 h=math.isqrt(budget//D[j]);off=sum(FF[j][i]*V[i] for i in range(j+1,7));D0=FF[j][j]
 for v in range(-((h+off)//D0),(h-off)//D0+1):
  V[j]=v;enum(j-1,budget-D[j]*(D0*v+off)**2)
enum(6,1440);assert np.max(tab)<=24
def costs(bs):
 return np.array([[tab[sum((((bb>>i&1)+2*(mm>>i&1))%4)<<(2*i) for i in range(7))] for mm in M] for bb in bs])
# An independent 4-block field solver: choose two field symbols and three radical bits; solve for other entries, with no min-plus convolution.
j=np.arange(32768);k0=j&63;k1=(j>>6)&63
h0=(j>>12)&1;h1=(j>>13)&1;h2=(j>>14)&1
inv=next(t for t in range(64) if fld(A[2]^A[3],t)==1)
def extmin(cost,s):
 a=(s&63)^k0^k1; b=(s>>6&63)^mul[A[0],k0]^mul[A[1],k1]
 k2=mul[inv,b^mul[A[3],a]];k3=a^k2;h3=(s>>12)^h0^h1^h2
 result=cost[0,k0+64*h0]+cost[1,k1+64*h1]+cost[2,k2+64*h2]+cost[3,k3+64*h3]
 return int(np.min(result))
def local(qq,m):
 # If squared norm <= 8, |2y_j-qq_j-m|<=5, so this box contains every vector needed.
 lower=[(q+m-6)//2 for q in qq]; upper=[(q+m+6)//2 for q in qq]; n=np.prod([h-l+1 for l,h in zip(lower,upper)])
 y=np.array(np.unravel_index(np.arange(n),[h-l+1 for l,h in zip(lower,upper)]),dtype=np.int16)+np.array(lower)[:,None]
 total=y.sum(axis=0);n2=((2*y-np.array(qq)[:,None]-m)**2).sum(axis=0);lab=(np.array([1,2,4,8,16])@(y[:5]%2))+32*((total//2)%2)
 return [min(32,int(n2[(total%2==0)&(lab==x)].min())) for x in range(64)]
def layers(q,b,u):
 out=[];hh=[]
 for m in range(2):
  bs=[z^(T[i] if m else 0) for i,z in enumerate(b)]
  # Coefficients b-r = (b XOR r)-2(r AND NOT b). Transfer the integral correction to the syndrome.
  si=u^(of if m else 0)
  if m:
   for i in range(4):si^=syn(i,T[i]&~b[i])
  co=costs(bs);lo=local(q,m); vals=[(lo[t],extmin(co,si^fx(t))) for t in range(64)]
  out.append(min(x+y for x,y in vals));hh.append(sorted(Counter(vals).items()))
 return out,hh
d1,hist=layers(hloc,bIn,u);print("first layer: four times squared distances to 34D parent banks",d1)
assert min(d1)==26
# Double the new shift. Its local coordinate hloc lies in D6. Remove that projected local fiber; doubling the chosen integral representative has zero syndrome.
p=sum((v%2)<<i for i,v in enumerate(hloc[:5]))+32*((sum(hloc)//2)%2)
ss=0
for i,b in enumerate(bIn):ss^=syn(i,b)
net=ss^fx(p); print("double: local tag, external tag, reduced syndrome",fx(p),ss,net);assert [fx(p),ss,net]==[88,919,975]
d2,hh=layers([0]*6,[0]*4,net);print("second layer: four times squared distances",d2);assert min(d2)>=8
print("Squared minimum 8, covolume 32768 sqrt(3/2), center density 8/sqrt(3).")
