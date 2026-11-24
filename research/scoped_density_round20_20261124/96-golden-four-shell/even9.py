from fractions import Fraction as F
from collections import defaultdict
from itertools import product
from pathlib import Path
import sympy as S
from anchored_replay import matrix
from jacobi_replay import J,K,E4,E6,D,powd,mul
from functools import lru_cache
@lru_cache(None)
def one(m,a):
 x=mul(powd(J,a),powd(K,m-a))
 return [{k:v for (i,k),v in x.items() if i==n} for n in range(9)]

def basis(m,l,a,b,j,w):
 assert w%2==0
 u=next((i,h) for i in range(w//4+1) for h in range(2) if 4*i+6*h==w)
 E=mul(powd(D,j),mul(powd(E4,u[0]),powd(E6,u[1])))
 A,B=one(m,a),one(l,b)
 out=[defaultdict(int) for _ in range(9)]
 for (i,z),v in E.items():
  assert not z
  for k in range(9-i):
   for h in range(9-i-k):
    for p,c in product(A[k], B[h]):
     if min(p,c)>=0:out[i+k+h][p,c]+=v*A[k][p]*B[h][c]
 return [dict(x) for x in out]

def staircase(m=5,l=5):
 L=list(product(range(m+1),range(l+1)))
 M=S.Matrix([[one(m,a)[0].get(p,0)*one(l,b)[0].get(c,0) for a,b in L] for p,c in L]).inv()
 out=[defaultdict(F) for i in range(9)]
 out[0][0,0]=-1
 for j in range(5):
  C=M*S.Matrix([-out[j][p,c] for p,c in L])
  for (a,b),t in zip(L,C):
   if not t: continue
   E=basis(m,l,a,b,j,48+2*(a+b)-12*j)
   for i in range(j,9):
    for k,v in E[i].items():out[i][k]+=F(t)*v
  assert not any(out[j].values())
 out[0][0,0]=1
 space=[basis(m,l,a,b,5,2*(a+b)-12) for a,b in L if a+b==6 or a+b>=8]
 return out,space

def evenrows():
 R,V,M,B,Fb=matrix()
 C,E=staircase()
 D=[(p,c) for p,c in product(range(11),repeat=2) if (p<=5 and c<=5) or (p,c)in((10,0),(0,10)) or C[5].get((p,c)) or any(e[5].get((p,c)) for e in E)]
 Q=S.Matrix([[e[5].get((p,c),0) for e in E]for p,c in D])
 print('even free',len(E),'rank',Q.rank())
 for v in Q.T.nullspace():
  vv=[F(z) for z in v]; hh={q:z for q,z in zip(D,vv) if z}
  def val(k):return hh.get(tuple(map(abs,k[:2])),F(0))/(2**sum(t!=0 for t in k[:2]))
  row=[sum(val(k) for k in bb) for bb in B]
  b=sum(t*C[5].get(k,0) for t,k in zip(vv,D))
  R.append(row);V.append(b);M.append(('even',len(R)))
 return R,V,M,B

if __name__=='__main__':
 import json
 R,V,M,B=evenrows()
 D={int(k):F(v) for k,v in json.loads(Path('anchor-counts.json').read_text()).items()}
 fails=[(m,sum(r[i]*v for i,v in D.items())-w) for r,w,m in zip(R,V,M) if sum(r[i]*v for i,v in D.items())!=w]
 print('old failures',len(fails),fails[:2])
