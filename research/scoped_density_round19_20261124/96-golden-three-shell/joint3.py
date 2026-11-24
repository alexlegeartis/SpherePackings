from second import *
from math import lcm,gcd

# Odd feature generator phi_{-1,2}, normalized with leading r-r^{-1}.
JJ={(0,1):1,(0,-1):-1}
for n in range(1,8):
 JJ=mul(JJ,mul(mul({(0,0):1,(n,2):-1},{(0,0):1,(n,-2):-1}),
                {(n*k,0):k+1 for k in range(8//n+1)}))
def odd1(a):
 out=mul(JJ,mul(powd(J,a),powd(K,3-a)))
 return [{p:w for (i,p),w in out.items() if i==n and p>0} for n in range(3)]
def odd():
 out=[]
 for a,b in ((2,3),(3,2)):
  H=[];A,B=odd1(a),odd1(b)
  for n in (0,1,2):
   C=defaultdict(int)
   for i in range(n+1):
    for j in range(n+1-i):
     W=(1,-120,7020)[n-i-j] # Delta^5
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
  B0=C[5:8] if not a else [{},{},{}]
  space=[e[5:8] for e in E] if not a else E2
  yield a,B0,space


from shell14 import orbits as orbits14, shellrows as rs14, N12 as N14
def jointrows3():
 from joint import jointrows
 R,V,M,A,B=jointrows()
 P,Fb=orbits14();assert not Fb
 FA,FB,FP=len(A),len(B),len(P)
 R=[r+[F(0)]*FP for r in R]
 H,W,T=rs14(P)
 R += [[F(0)]*(FA+FB)+r for r in H];V+=W;M+=T
 for od,C,E in spaces():
  K=range(1 if od else 0,12)
  D=sorted(set(product(K,repeat=2)))
  Q=S.Matrix([[e[0].get(k,0) for e in E]for k in D])
  sel=Q.T.rref()[1]; QQ=Q[list(sel),:]
  assert len(sel)==len(E)
  PR=[parity(A,D[i],od) for i in sel]
  PV=[C[0].get(D[i],0) for i in sel]
  QI=QQ.inv()
  for k in D:
   rr=[F(x) for x in S.Matrix([[e[2].get(k,0) for e in E]])*QI]
   wr=parity(P,k,od)
   if not any(rr) and not any(wr):
    assert not C[2].get(k,0);continue
   R.append([-sum(t*r[i] for t,r in zip(rr,PR))for i in range(FA)]+[0]*FB+wr)
   V.append(C[2].get(k,0)-sum(t*v for t,v in zip(rr,PV)))
   M.append(('joint theta14',od,k))
 lookup={q:i for i,b in enumerate(P)for q in b}
 for X,d,f in ((A,10,0),(B,12,FA)):
  for i,bb in enumerate(X):
   for p,c,s in bb:
    opts=[]
    if (3*d+s)//2==14:opts.append((c,p+c,(5*d+3*s)//2))
    if (3*d-s)//2==14:opts.append((c-p,p,(-5*d+3*s)//2))
    for u,v in product(range(-1,2),repeat=2):
     if not(u or v) or 5*(u*u+v*v)+p*u+c*v!=(14-d)//2:continue
     h=-10*u*u+40*u*v+10*v*v+2*((2*c-p)*u+(2*p+c)*v)
     opts.append((p+10*u,c+10*v,s+h))
    for z in opts:
     out=[0]*(FA+FB+FP);out[f+i]=1
     if z in lookup:out[FA+FB+lookup[z]]=-1
     R.append(out);V.append(0);M.append(('cross14',f+i,z))
 return R,V,M,A,B,P
if __name__=='__main__':
 rr=jointrows3()
 print(len(rr[0]),[len(t)for t in rr[3:]])
