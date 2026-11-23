from pathlib import Path
from itertools import product
import json,sys,subprocess,tempfile
import numpy as np
from sympy import Matrix
P=Path(__file__).resolve().parent
C=json.loads((P/'certificate.json').read_text())
M=C['phase_isometry']; alpha=C['alpha']
def mul(a,b):
 r=0
 while b:
  if b&1:r^=a
  b>>=1;a<<=1
  if a&64:a^=67
 return r
def rep(s):
 x,y,e=s&63,s>>6&63,s>>12
 return [M[(x^y)|(e<<6)],M[y],0,0,0]
def unpack(a):return [(a>>(2*i))&3 for i in range(7)]
def bits(a):return [(a>>i)&1 for i in range(7)]
R=bits(82)
B=[272,1360,1,1089,0]
S=[439,291,684]
A2=np.array([[R+np.zeros(7,dtype=int)+2*np.array(bits(x)) for x in rep(S[0])],
             [np.array(unpack(a))+2*np.array(bits(x))for a,x in zip(B,rep(S[1]))],
             [np.array(unpack(a))+2*np.array(bits(x))for a,x in zip(B,rep(S[2]))]],dtype=int)
V2=np.array([[1,1,0],[2,0,2],[0,2,2]],dtype=int) # twice vertical
assert Matrix(V2.tolist()).det()==-8 or Matrix(V2.tolist()).det()==8
inv={v:i for i,v in enumerate(M)}
def syn(v):
 s=0
 for i,o in enumerate(v):
  y=inv[sum((int(w)&1)<<j for j,w in enumerate(o))]
  s^=(y&63)|(mul(alpha[i],y&63)<<6)|((y>>6)<<12)
 return s
assert syn(2*A2[0])==0
assert syn(A2[1])==0 and syn(A2[2])==0
# 4 times first horizontal shift and twice either other shift lie in L35.
G=Matrix(C['E7_gram']); L,D=G.LDLdecomposition(hermitian=False)
raw=list(map(int,(P/'local4.in').read_text().split()));s0,t0=raw[:2]
assert all(s0*L[i,j]==raw[9+7*i+j]for i in range(7)for j in range(7)) and all(t0*D[i,i]==raw[2+i]for i in range(7))
with tempfile.TemporaryDirectory() as t:
 b=Path(t)/'local4'
 subprocess.run(['g++','-O3','-std=c++17',str(P/'local4.cpp'),'-o',str(b)],check=True)
 z=subprocess.run([str(b)],input=(P/'local4.in').read_bytes(),capture_output=True,check=True)
D4=np.fromstring(z.stdout.decode(),sep=' ',dtype=np.int32)
assert len(D4)==16384 and max(D4)==24
cols=[[int((x&63)|(mul(alpha[j],x&63)<<6)|((x>>6)<<12))for x in range(128)]for j in range(5)]
Q=[sum(int(w)<<(2*j)for j,w in enumerate(bits(x)))for x in M]
ids=np.arange(8192,dtype=np.int32)
ix=[ids[None,:]^np.asarray(row,dtype=np.int32)[:,None]for row in cols]
D=[]
for k,l,m in product(range(4),range(2),range(2)):
 z=(np.array([k,l,m],dtype=int)@V2).tolist()
 vv=min(sum((a+4*b)**2 for a,b in zip(z,n))for n in product(range(-4,5),repeat=3)if sum(n)%2==0)
 assert vv<=32 # any omitted n_i has |z_i+4n_i|>sqrt(32)
 X=k*A2[0]+l*A2[1]+m*A2[2]
 f=np.full(8192,100,dtype=np.int32);f[0]=0
 for j,y in enumerate(X):
  raw=sum((int(a)&3)<<(2*i)for i,a in enumerate(y))
  u=D4[np.asarray([raw^(2*p)for p in Q],dtype=np.int32)]
  f=(u[:,None]+f[ix[j]]).min(axis=0)
 h=int(f[0])
 if k or l or m:assert h+vv>=32, (k,l,m,h,vv)
 D.append([k,l,m,h,vv])
Cout=dict(beta_packed=B,syndrome_shifts=S,horizontal_generators_times_2=A2.tolist(),
          vertical_generators_times_2=V2.tolist(),orbit_costs_times_4=D,
          gram_determinant=2**31,squared_minimum=8)
if '--write' in sys.argv:(P/'certificate38.json').write_text(json.dumps(Cout,indent=2)+'\n')
else:assert Cout==json.loads((P/'certificate38.json').read_text())
for t in D:print(*t)
print('PASS exact 16 horizontal/vertical residue costs, determinant 2^31, min squared 8, center density 8 sqrt(2)')
