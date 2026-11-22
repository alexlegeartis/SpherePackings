import re,json,itertools,hashlib,math,sys
from pathlib import Path
from fractions import Fraction as F
from flint import fmpz_mat as Z
P=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent.parent/'packing78_verified_20261117'/'packing78_independent_20261117'/'L80_source.html'
b=P.read_bytes()
assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()=='ebab4434e56f176cbe2ae6c8346a90b1cd05245b'
s=b.decode().split('NAME="GRAM"')[1].split('NAME="DET"')[0]
n,k,*z=map(int,re.findall(r'-?\d+',re.sub('<[^>]*>',' ',s))); assert n==80 and k==0 and len(z)==3240
G=[[0]*n for _ in range(n)]; t=0
for i in range(n):
 for j in range(i+1):G[i][j]=G[j][i]=z[t];t+=1
G=Z(G); assert G.det()==1
c=[i for i in range(n)if G[i,i]==8]
for i,j,k in itertools.permutations(c,3):
 if abs(G[i,j])==abs(G[j,k])==4 and G[i,k]==0:break
else:raise AssertionError('no short catalogue triple')
s=[-int(G[i,j])/4,1,-int(G[j,k])/4]
s=list(map(int,s)); assert set(s)<={-1,1}
R=Z([[sign*int(x==a)for x in range(n)]for a,sign in zip([i,j,k],s)])
A=R*G*R.transpose();assert A.tolist()==[[8,-4,0],[-4,8,-4],[0,-4,8]] and A.det()==256
W=(R*G).transpose()
H,V=W.hnf(transform=True);assert abs(V.det())==1 and H.tolist()==[[int(i==j)for j in range(3)]for i in range(80)]
U=Z(V.tolist()[3:]); assert (U*W).is_zero() and (U*G*U.transpose()).det()==256
Z0=None
for i0 in range(80):
 for j0 in range(i0+1,80):
  for eps in (-1,1):
   if G[i0,i0]+G[j0,j0]+2*eps*G[i0,j0]==8 and all(W[i0,k]+eps*W[j0,k]==0 for k in range(3)):
    Z0=[int(i==i0)+eps*int(i==j0)for i in range(80)];break
  if Z0 is not None:break
 if Z0 is not None:break
if Z0 is None:
 for i0,j0,k0 in itertools.combinations(range(80),3):
  for e,f in itertools.product((-1,1),repeat=2):
   if G[i0,i0]+G[j0,j0]+G[k0,k0]+2*(e*G[i0,j0]+f*G[i0,k0]+e*f*G[j0,k0])==8 and all(W[i0,k]+e*W[j0,k]+f*W[k0,k]==0 for k in range(3)):
    Z0=[int(i==i0)+e*int(i==j0)+f*int(i==k0)for i in range(80)];break
  if Z0 is not None:break
assert Z0 is not None
VV=Z([Z0]);assert (VV*W).is_zero() and (VV*G*VV.transpose())[0,0]==8
def p(k,n):return sum(F((-1)**j,(2*j+1)*k**(2*j+1))for j in range(n+1))
lo=16*p(5,71)-4*p(239,70)
vol2=F(2**78,math.prod(range(1,78,2))**2)*lo**76
dn=F('1.03072853214953e-16')
# Print the conservative floor corresponding to the exact squared center density.
from decimal import Decimal,localcontext,ROUND_FLOOR
with localcontext() as ctx:
 ctx.prec=100
 sq=F(2**69)*vol2
 d=(Decimal(sq.numerator)/Decimal(sq.denominator)).sqrt().quantize(Decimal('1e-30'),rounding=ROUND_FLOOR)
assert F(d)**2<sq
print(i,j,k,s,Z0,'det256 minimum8 exact; occupied >=',d)
Path('published77.json').write_text(json.dumps({'basis_indices_zero_based':[i,j,k],'signs':s,'minimum8_witness':Z0,'normal_gram':A.tolist(),'determinant':256,'center_squared':'2^69','occupied_floor':str(d)},default=int,indent=2)+'\n')
