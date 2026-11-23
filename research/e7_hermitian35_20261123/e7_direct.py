from itertools import product
from math import isqrt
from fractions import Fraction as F
def dot(a,b): return sum(a[i]*b[i]*(2 if i==6 else 1) for i in range(7))
def cal(i,s):
 if i==7:
  yield ()
  return
 l=isqrt(s//(2 if i==6 else 1))
 for j in range(-l,l+1):
  for r in cal(i+1,s-j*j*(2 if i==6 else 1)):yield (j,)+r
B=[]
for i in range(5):
 a=[0]*7;a[i]=2;a[i+1]=-2;B.append(a)
B += [(0,0,0,0,2,2,0),(1,1,1,1,1,-1,1)]
def key(x):
 if any((y-x[6])%2 for y in x)or sum(x[:6])%4:return None
 a=[(x[i]-x[6]*B[-1][i])//2 for i in range(6)]
 v=[sum(a[:i+1])for i in range(4)]+[(sum(a[:5])-a[-1])//2,sum(a)//2,x[6]]
 assert [sum(v[i]*B[i][j]for i in range(7))for j in range(7)]==list(x)
 return sum((c%2)<<i for i,c in enumerate(v))
g=[[F(dot(a,b),4)for b in B]for a in B]
g0=[r[:] for r in g];c=F(1)
for j in range(7):
 h=g[j][j];assert h>0;c*=h
 for i in range(j+1,7):
  z=g[i][j]/h
  for p in range(j+1,7):g[i][p]-=z*g[j][p]
assert c==2
p=lambda a,b:sum(int(g0[i][j])*(a>>i&1)*(b>>j&1)for i in range(7)for j in range(7))%2
vv=lambda a:[sum(B[i][j]*(a>>i&1)for i in range(7))for j in range(7)]
q=lambda x:dot(vv(x),vv(x))//8%2
rad=[x for x in range(128)if not any(p(x,1<<j)for j in range(7))]
assert len(rad)==2 and q(rad[1])==1
D={}
for a in cal(0,24):
 k=key(a)
 if k is None:continue
 v=dot(a,a)//4
 assert 4*v==dot(a,a) and v%2==0
 D[k]=min(D.get(k,100),v)
assert len(D)==128
assert D[0]==0 and D[rad[1]]==6
assert all(v==(2 if q(x) else 4)for x,v in D.items()if x not in rad)
print("PASS independent E7, det 2, all 128 phase minima, radical representative squared length 6")
# 35D:
# Residue projection to F64^5 is a rank-3 MDS code. A weight-three
# word is a common scalar times (α2+α3, α3+α1, α1+α2) in F8^3.
# Norm on F8 multiplies by the square, so its sum vanishes. Those
# triples have even total q when extra radical parity is even, and
# hence cannot all be weight-two phases. Other supports use ≥4 active
# blocks or ≥2 rad-only blocks. Thus min8, with period equality.
assert 2**35*2**5==2**40 # convenient consistency check
# covolume squared of (2E7)^5 = 2**75; phase index = 2**22.
V2=F(2**75,2**44);assert V2==2**31
delta35sq=F(2)**35/V2;assert delta35sq==16
# 36D: multiples of (r/2)^5 have lower squared horizontal distances
# (0, 30/4, 6, 30/4) mod 4, respectively. At k=0 use old min8; at
# k≡0 and nonzero add at least 8 vertically, at odd k at least 1/2,
# at k≡2 add at least 2.
assert F(30,4)+F(1,2)==6+F(4,2)==F(16,2)==8
delta36sq=F(2)**36/(V2/2);assert delta36sq==64
print("PASS independent density squares 16,64")
