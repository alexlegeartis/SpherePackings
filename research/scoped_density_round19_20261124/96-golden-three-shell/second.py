from even import *
from anchored_replay import N
import json
from math import prod
from anchored_replay import COEF,q, SZ as Z10
S12={0:1218201435360000,4:726542605152000,8:143169114816000,12:7191815904000,16:49886928000}
S12|={-k:v for k,v in S12.items()}
N12=sum(S12.values())
def valid(p,c,s):
 m=(12+s)//2; n=(36+s)//2
 G=((10,0,p,c),(0,10,c,p+c),(p,c,12,m),(c,p+c,m,n))
 a=120-p*p-c*c;b=10*n-c*c-(p+c)**2;d=10*m-p*c-c*(p+c)
 return a>=0 and b>=0 and a*b>=d*d and not any(0<q(G,t)<10 for t in COEF)
def images(p,c,s):
 out=[(-p,-c,s)]
 if s==-12:out.append((c,p+c,12))
 if s==12:out.append((c-p,p,-12))
 for u,v in product(range(-1,2),repeat=2):
  if not(u or v) or 5*(u*u+v*v)+p*u+c*v: continue
  h=-10*u*u+40*u*v+10*v*v+2*((2*c-p)*u+(2*p+c)*v)
  out.append((p+10*u,c+10*v,s+h))
 return out
def orbits():
 A={(p,c,s)for p,c,s in product(range(-6,7),range(-6,7),S12) if valid(p,c,s)}
 B=[];forbid=[]
 while A:
  x=min(A);todo=[x];Z={x}
  while todo:
   for y in images(*todo.pop()):
    if y not in Z:todo.append(y);Z.add(y)
  (B if Z<=A else forbid).append(tuple(sorted(Z&A)))
  A-=Z
 return B,forbid

def shellrows(B):
 R=[];V=[];M=[]
 def add(f,t,l):
  R.append([sum(f(*a) for a in bb) for bb in B]);V.append(F(t));M.append(l)
 for s,A in S12.items():
  add(lambda p,c,h:int(h==s), A,('size',s))
  add(lambda p,c,h:(5*p*p+(2*c-p)**2)*int(h==s),F(600-10*s,48)*A,('cov1',s))
  add(lambda p,c,h:p*(2*c-p)*int(h==s),F(5*s-60,48)*A,('cov2',s))
 for i in range(11):
  for j in range(11-i):
   if (i+j)%2 or not (i+j):continue
   t=0
   if not i%2 and not j%2:
    m=(i+j)//2
    t=F(N12*120**m*prod(range(1,i,2))*prod(range(1,j,2)),prod(96+2*h for h in range(m)))
   add(lambda p,c,s:p**i*c**j,t,('mom',i,j))
 return R,V,M
def predicted(C0={}):
 R,V,M,B=evenrows()
 if not C0: C0={int(k):F(v) for k,v in json.loads(Path('anchor-counts.json').read_text()).items()}
 C,E=staircase()
 D=[(p,c) for p,c in product(range(12),repeat=2) if (p<=6 and c<=6) or (p,c)in((10,0),(0,10)) or C[5].get((p,c)) or C[6].get((p,c)) or any(e[5].get((p,c)) or e[6].get((p,c)) for e in E)]
 def vals(k,B):
  return sum(F(sum(tuple(map(abs,z[:2]))==k for z in bb),2**sum(t!=0 for t in k))*C0.get(j,0) for j,bb in enumerate(B))
 Q=S.Matrix([[e[5].get((p,c),0) for e in E]for p,c in D])
 x=Q.solve(S.Matrix([vals((p,c),B)-C[5].get((p,c),0) for p,c in D]))
 return {k:C[6].get(k,0)+sum(e[6].get(k,0)*F(z)for e,z in zip(E,x)) for k in D}
if __name__=='__main__':
 Z=predicted()
 print('q6 neg',[(k,v)for k,v in Z.items()if v<0])
 print('q6 beyond',[(k,v)for k,v in Z.items()if max(k)>6 and v])
 B,Fb=orbits();print('orbits',len(B),'invalid',len(Fb))
 print('q6 support valid',all(any(tuple(map(abs,y[:2]))==k for bb in B for y in bb) or not v for k,v in Z.items()))
 print('q6', [(k,float(v)) for k,v in Z.items() if v and max(k)>=5])
