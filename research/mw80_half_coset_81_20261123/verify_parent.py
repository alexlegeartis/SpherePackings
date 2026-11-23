from pathlib import Path
import json,numpy as np
from flint import fmpz_mat as Z,fq_default_ctx,fq_default_poly_ctx,fmpz_mod_poly_ctx
p=41;F=fq_default_ctx(modulus=fmpz_mod_poly_ctx(p)([-3,0,1]),fq_type='FQ_NMOD')
R=fq_default_poly_ctx(F); f=R([1]+[0]*41+[1])
A=np.load('geometry.npz');Y=A['roots'];G=A['G'];ZZ=Z(G.tolist());YY=Z(Y.tolist())
assert int(ZZ.det())==p**12
D=json.loads(Path('unitaries80.json').read_text());assert D['prime']==41 and D['field_nonsquare']==3
assert len(D['unitaries'])==len(Y)
xs=[];ys=[];keys=set()
for row in D['unitaries']:
 a,b,c,d=[F([z%p,z//p])for z in row]
 assert a**42+b**42==c**42+d**42==1 and a*c**41+b*d**41==0
 x=-(R([b,a])**14);y=R([d,c])**21
 assert y*y==x**3+f and (str(x),str(y))not in keys
 keys.add((str(x),str(y)));xs.append(x);ys.append(y)
H=YY*ZZ*YY.transpose();n=len(xs)
for i,(x,y) in enumerate(zip(xs,ys)):
 assert H[i,i]==14
 for j in range(i):
  xx,yy=x-xs[j],y-ys[j]
  if xx.is_zero() and (y+ys[j]).is_zero():v=-14
  else:
   a=xx.gcd(yy).degree()
   b=min(([14-xx.degree()]if not xx.is_zero()else[])+([21-yy.degree()]if not yy.is_zero()else[]))
   v=7-a-b
  assert H[i,j]==v
assert YY.hnf().tolist()[:80]==Z(np.eye(80,dtype='i8').tolist()).tolist()
print('Exact p41 MW80 height matrix, primitive full row-generated basis, det',p**12,'min14 via Shioda height formula')
