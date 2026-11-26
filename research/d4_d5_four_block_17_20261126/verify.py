"""Standalone exact 17D four-block certificate. No dependencies, network or optimizer."""
from fractions import Fraction as F
from functools import reduce
from itertools import product, permutations, combinations, combinations_with_replacement as cwr
from math import prod,comb
from pathlib import Path
from json import loads
from operator import xor
D=Path(__file__).resolve().parent
H=[[1,1,1],[1,-F(1,3),1],[1,0,-F(1,3)]]
O=[[1,1,1,1],[1,F(1,5),-F(3,5),1],[1,-F(1,5),F(1,5),1],[1,0,0,-1]]
# Independent Fourier tables from binary characters (D4* alphabet and D5).
def chi(a,x):return (-1)**((a&x).bit_count())
def qc(x):return 0 if x==0 else 2 if x<4 else 1
def oc(x):
 p=x%16;w=p.bit_count();w+=w%2
 return 3 if not p and x else w//2
def ft(n,types,cl):
 sets=[[x for x in range(n) if cl(x)==t]for t in range(types)]
 return {tuple(F(sum(chi(a,x) for x in z),len(z))for z in sets)for a in range(n)}
assert ft(16,3,qc)==set(map(tuple,H)) and ft(32,4,oc)==set(map(tuple,O))
def check_orbit_means(n,cl,R,m):
 groups=[[] for _ in R];sets=[[x for x in range(n)if cl(x)==t]for t in range(len(R[0]))]
 for a in range(n):
  row=[F(sum(chi(a,x)for x in z),len(z))for z in sets]
  groups[R.index(row)].append(a)
 assert [len(orb)for orb in groups]==m
 for row,orb in zip(R,groups):
  assert all(F(sum(chi(a,x)for a in orb),len(orb))==row[cl(x)]for x in range(n))
check_orbit_means(16,qc,H,[1,3,12]);check_orbit_means(32,oc,O,[1,5,10,16])
ids=list(product(range(3),range(3),range(3),range(4)))
prof=[(a,b,c,[0,1,2,2][d]) for a,b,c,d in ids]
par=[(p,4) for p in sorted(set(prof)) if 1<=sum(p)<=3]+[(tuple(int(i==j) for i in range(4)),1) for j in range(4)]
def det3(a,b,c):return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def det(M):return sum((-1)**i*M[0][i]*det3(*[a[:i]+a[i+1:] for a in M[1:]]) for i in range(4) if M[0][i])
def vertices():
 out=set()
 for pp in combinations(par,4):
  M=tuple(p for p,b in pp);b=tuple(b for p,b in pp); d=det(M)
  if not d:continue
  r=tuple(F(det(tuple(row[:i]+(bi,)+row[i+1:] for row,bi in zip(M,b))),d) for i in range(4))
  if min(r)>=1 and max(r)<=4:out.add(tuple(sorted(r[:3]))+r[3:])
 return out
K=[[H[x[0]][y[0]]*H[x[1]][y[1]]*H[x[2]][y[2]]*O[x[3]][y[3]] for y in ids]for x in ids]
W=loads((D/'scales.json').read_text());W={tuple(map(F,k.split(','))):v for k,v in W.items()}
assert set(W)==vertices() and len(W)==209
for r,ww in W.items():
 assert all(0<i<len(ids) and w>=0 for i,w in ww)
 S=10**6;I=[j for j,p in enumerate(prof) if sum(x*y for x,y in zip(p,r))>=4]
 gap=min(-sum(w*K[i][j] for i,w in ww) for j in I)
 assert gap>0
 U=1+F(sum(w for i,w in ww),gap)
 if r==(1,1,1,1):
  assert int(U)==397
 else:
  assert int(U)**2 <= 2*256**2*prod(r[:3])**4*r[3]**5
# Auxiliary even-word lift: all distances and all characters, with no affine-independence assumption.
J=list(product(cwr(range(3),3),range(4)));u=[0,2,4,6]
A=[[sum(prod(H[i][k] for i,k in zip(p,y[0])) for p in permutations(x[0]))*O[x[1]][y[1]]/F(6) for y in J]for x in J]
B=[[F(sum((-1)**j*comb(s,j)*comb(6-s,v-j)for j in range(min(v,s)+1)if 0<=v-j<=6-s),comb(6,v)) for v in u]for s in range(4)]
du=loads((D/'dual6.json').read_text());y,z,S=du['y'],du['z'],du['quant'];assert len(y)==160 and all(v>=0 for v in y) and len(z)==3
eps=F(0)
for j,x in enumerate(J):
 for iv,v in enumerate(u):
  if v<=2 and sum(x[0])+[0,1,2,2][x[1]]<4:continue
  excess=sum(Y*A[i//4][j]*B[i%4][iv]for i,Y in enumerate(y)if Y) - sum(z[ik-1]*(int(iv==ik)-comb(6,u[ik])*int(iv==0))for ik in range(1,4))+S*int(iv==0)
  eps=max(eps,excess)
q=[15,15,1];assert S>32*eps
bound=F(S+sum(y)+sum(aa*bb for aa,bb in zip(q,z))-eps,S-32*eps)
assert bound==F(100517603567,1874999344) and bound<54
gs=[31181,58116,71012,35581,40360,35279,18,258]
words={0}
for v in gs:words|={x^v for x in words}
def dist(w):return sum(qc((w>>k)&15)for k in (0,4,8))+[0,1,2,2][oc(w>>12)]
assert len(words)==256 and min(map(dist,words-{0}))==4
print("Verified 209 shapes, 208 excluded; isotropic M <=397; six-packet kernel <=",bound,"; positive 256 centers, not a record.")

