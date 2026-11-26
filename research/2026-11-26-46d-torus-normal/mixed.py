from fractions import Fraction as F
from math import factorial,prod,lcm,gcd
from itertools import product
from functools import cache
from collections import defaultdict
from sympy.polys.matrices import DomainMatrix
from sympy import ZZ
from sys import argv
import json
def exact(A,B,g,Dmax=16):
 D=A*B-g*g
 def key(a,b):
  k=((B*a-g*b)%D,(A*b-g*a)%D)
  return min(k,((-k[0])%D,(-k[1])%D))
 K=defaultdict(list)
 for a,b in product(range(-12,13),repeat=2):
  q=F(B*a*a-2*g*a*b+A*b*b,D)
  if q<=12:K[key(a,b)].append((q,a,b))
 keys=sorted(K); n=len(keys)
 tr=[]
 for k in keys:
  s=min(q for q,a,b in K[k])
  tr.append([[ (a,b) for q,a,b in K[k] if q==s+j] for j in (0,2)])
 pp={m:[(A*i+g*j,g*i+B*j) for i,j in product(range(-3,4),repeat=2) if A*i*i+2*g*i*j+B*j*j==m] for m in (6,8)}
 @cache
 def harmonic(i,j,m):
  d=i+j
  P={(i,j):1}; U=defaultdict(F);c=F(1)
  for h in range(d//2+1):
   Q=defaultdict(int)
   for (a,b),v in P.items():
    U[a,b]+=c*v
    if a>=2:Q[a-2,b]+=v*a*(a-1)*A
    if b>=2:Q[a,b-2]+=v*b*(b-1)*B
    if a*b:Q[a-1,b-1]+=v*a*b*2*g
   P=Q;c*=F(-m,2*(h+1)*(2*d+48-2*(h+1)-2))
  t=lcm(*(v.denominator for v in U.values()))
  return t,[(a,b,int(v*t)) for (a,b),v in U.items() if v]
 def row(i,j,m):
  t,P=harmonic(i,j,m)
  @cache
  def f(x,y):return sum(v*x**a*y**b for a,b,v in P)
  get=lambda marks:sum(f(a,b) for a,b in marks)
  if m==6:R=[get(r[0]) for r in tr]+[0]*n
  else:R=[get(r[1]) for r in tr]+[get(r[0]) for r in tr]
  R+=[(0 if i+j else {6:52416000,8:39007332000}[m])*t-get(pp[m])]
  return R,t
 eq=[]
 for d in range(0,Dmax+1,2):
  for i in range(d+1):
   j=d-i
   x,t=row(i,j,6);y,u=row(i,j,8)
   if d<=10 or d==14:RR=[x,y]
   else:
    lam={12:-72,16:168,18:-576,20:408,22:-336}[d]
    RR=[[b*t-lam*a*u for a,b in zip(x,y)]]
   for r in RR:
    h=gcd(*r)
    if h:eq.append([v//h for v in r])
 M=DomainMatrix.from_list(eq,ZZ);M,den,pivs=M.rref_den()
 assert 2*n not in pivs
 free=[i for i in range(2*n) if i not in pivs]
 rem={}; rows=M.to_list(); fs=set(free)
 for k,r in zip(pivs,rows):
  rem[k]=[F(int(r[-1]),int(den))]+[F(-int(r[f]),int(den)) for f in free]
 for i,f in enumerate(free):rem[f]=[F(0)]+[int(i==j) for j in range(len(free))]
 print('pair',A,B,g,'classes',n,'eq/rank/free',len(eq),len(pivs),len(free),flush=True)
 out={'A':A,'B':B,'g':g,'free':free,'keys':keys,'minima':[str(min(q for q,a,b in K[k])) for k in keys],'coords':[r[0] for r in tr],'aff':[[str(x) for x in rem[i]] for i in range(2*n)]}
 with open(f'mixed_{A}_{B}_{g}.json','w') as f:json.dump(out,f)
if __name__=='__main__':
 for g in map(int,argv[3:]):exact(int(argv[1]),int(argv[2]),g)
