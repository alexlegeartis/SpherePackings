from second import *
from math import lcm,gcd

# Odd feature generator phi_{-1,2}, normalized with leading r-r^{-1}.
JJ={(0,1):1,(0,-1):-1}
for n in range(1,8):
 JJ=mul(JJ,mul(mul({(0,0):1,(n,2):-1},{(0,0):1,(n,-2):-1}),
                {(n*k,0):k+1 for k in range(8//n+1)}))
def odd1(a):
 out=mul(JJ,mul(powd(J,a),powd(K,3-a)))
 return [{p:w for (i,p),w in out.items() if i==n and p>0} for n in range(2)]
def odd():
 out=[]
 for a,b in ((2,3),(3,2)):
  H=[];A,B=odd1(a),odd1(b)
  for n in (0,1):
   C=defaultdict(int)
   for i in range(n+1):
    for j in range(n+1-i):
     W=1 if n-i-j==0 else -120 # Delta^5
     for p,c in product(A[i],B[j]):C[p,c]+=W*A[i][p]*B[j][c]
   H.append(dict(C))
  out.append(H)
 return out
def parity(B,k,odd=False):
 out=[]
 for bb in B:
  s=0
  for p,c,t in bb:
   if (abs(p),abs(c))==k:s+=(-1 if p*c<0 else 1) if odd else 1
  out.append(F(s,2**sum(x!=0 for x in k)))
 return out

def spaces():
 C,E=staircase()
 E2=odd()
 for a in (False,True):
  B0=C[5:7] if not a else [{},{}]
  space=[e[5:7] for e in E] if not a else E2
  yield a,B0,space

def jointrows():
 R,V,M,A=evenrows()
 B,Fb=orbits();assert not Fb
 FA,FB=len(A),len(B)
 R=[r+[F(0)]*FB for r in R]
 H,W,T=shellrows(B)
 R += [[F(0)]*FA+r for r in H];V+=W;M+=T
 for od,C,E in spaces():
  P=range(1 if od else 0,11)
  D=sorted(set(product(P,repeat=2)))
  Q=S.Matrix([[e[0].get(k,0) for e in E]for k in D])
  sel=Q.T.rref()[1]; QQ=Q[list(sel),:]
  assert len(sel)==len(E)
  PR=[parity(A,D[i],od) for i in sel]
  PV=[C[0].get(D[i],0) for i in sel]
  QI=QQ.inv()
  for k in D:
   rr=[F(x) for x in S.Matrix([[e[1].get(k,0) for e in E]])*QI]
   wr=parity(B,k,od)
   if not any(rr) and not any(wr):
    assert not C[1].get(k,0);continue
   R.append([-sum(t*r[i] for t,r in zip(rr,PR))for i in range(FA)]+wr)
   V.append(C[1].get(k,0)-sum(t*v for t,v in zip(rr,PV)))
   M.append(('joint theta',od,k))
 lookup={q:i for i,b in enumerate(B)for q in b}
 for i,bb in enumerate(A):
  for p,c,s in bb:
   opts=[]
   if s==-6:opts.append((c,p+c,16))
   if s==6:opts.append((c-p,p,-16))
   for u,v in product(range(-1,2),repeat=2):
    if not(u or v) or 5*(u*u+v*v)+p*u+c*v!=1:continue
    h=-10*u*u+40*u*v+10*v*v+2*((2*c-p)*u+(2*p+c)*v)
    opts.append((p+10*u,c+10*v,s+h))
   for z in opts:
    out=[0]*(FA+FB);out[i]=1
    if z in lookup:out[FA+lookup[z]]=-1
    R.append(out);V.append(0);M.append(('cross',i,z))
 return R,V,M,A,B

if __name__=='__main__':
 R,V,M,A,B=jointrows()
 print('rows',len(R),'cols',len(A),len(B))
