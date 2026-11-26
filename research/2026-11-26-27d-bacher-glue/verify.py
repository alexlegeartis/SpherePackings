"""Exact one-core replay of the 22+4 splitting and 27D phase costs. Python 3.10+; standard library only."""
from fractions import Fraction as F
from math import lcm,isqrt
from itertools import product
from collections import Counter
G27=[[8,4,4,4,-2,2,-4,-2,-4,-4,-4,0,-2,2,2,0,2,-4,-4,-2,-4,2,4,4,-2,2,4],[4,8,4,4,-4,-2,-4,2,-2,-4,-4,0,-2,-2,0,2,4,-4,0,-4,0,-2,2,4,2,-2,2],[4,4,8,4,-4,-2,-2,2,-2,-2,-2,-2,-2,-2,-2,-2,2,-4,0,-2,-2,2,0,4,-2,2,0],[4,4,4,8,0,2,0,0,-4,-2,-4,2,-4,-2,-2,0,0,-2,0,0,0,0,0,4,-2,0,0],[-2,-4,-4,0,8,4,4,-4,2,4,0,0,-2,2,-2,-2,-4,2,-2,4,0,-2,-2,-4,-2,2,0],[2,-2,-2,2,4,8,0,-2,0,2,0,0,0,4,0,-2,-4,2,-2,4,0,2,2,0,-4,0,1],[-4,-4,-2,0,4,0,8,-2,0,4,2,2,0,-2,-4,-2,-2,2,0,4,0,0,-4,-4,-2,2,-2],[-2,2,2,0,-4,-2,-2,8,0,0,0,0,2,-2,0,0,0,2,4,-2,2,0,-2,2,2,-4,-4],[-4,-2,-2,-4,2,0,0,0,8,4,4,-2,0,2,-2,0,0,2,2,0,2,0,-2,-4,0,0,-1],[-4,-4,-2,-2,4,2,4,0,4,8,4,-2,2,2,-4,-4,-4,4,2,4,2,0,-2,-4,-2,2,-3],[-4,-4,-2,-4,0,0,2,0,4,4,8,-2,2,0,0,-2,-2,4,2,2,2,2,0,-2,0,0,-3],[0,0,-2,2,0,0,2,0,-2,-2,-2,8,-2,0,0,4,2,2,0,-2,0,0,-2,0,0,-2,-1],[-2,-2,-2,-4,-2,0,0,2,0,2,2,-2,8,2,2,0,-2,2,2,2,0,2,2,-2,0,0,-1],[2,-2,-2,-2,2,4,-2,-2,2,2,0,0,2,8,2,0,-2,2,-2,0,-2,2,2,-2,-2,2,2],[2,0,-2,-2,-2,0,-4,0,-2,-4,0,0,2,2,8,2,0,0,-2,-2,-2,0,4,2,2,-2,1],[0,2,-2,0,-2,-2,-2,0,0,-4,-2,4,0,0,2,8,4,0,2,-4,0,0,0,0,2,-2,1],[2,4,2,0,-4,-4,-2,0,0,-4,-2,2,-2,-2,0,4,8,-4,0,-4,-2,0,0,2,2,0,3],[-4,-4,-4,-2,2,2,2,2,2,4,4,2,2,2,0,0,-4,8,2,2,2,0,-2,-2,0,-2,-5],[-4,0,0,0,-2,-2,0,4,2,2,2,0,2,-2,-2,2,0,2,8,0,4,0,-2,0,2,-2,-3],[-2,-4,-2,0,4,4,4,-2,0,4,2,-2,2,0,-2,-4,-4,2,0,8,0,0,0,-2,-2,2,-1],[-4,0,-2,0,0,0,0,2,2,2,2,0,0,-2,-2,0,-2,2,4,0,8,-2,0,0,2,-4,-2],[2,-2,2,0,-2,2,0,0,0,0,2,0,2,2,0,0,0,0,0,0,-2,8,0,0,-4,2,1],[4,2,0,0,-2,2,-4,-2,-2,-2,0,-2,2,2,4,0,0,-2,-2,0,0,0,8,2,0,0,3],[4,4,4,4,-4,0,-4,2,-4,-4,-2,0,-2,-2,2,0,2,-2,0,-2,0,0,2,8,0,-2,0],[-2,2,-2,-2,-2,-4,-2,2,0,-2,0,0,0,-2,2,2,2,0,2,-2,2,-4,0,0,8,-2,0],[2,-2,2,0,2,0,2,-4,0,2,0,-2,0,2,-2,-2,0,-2,-2,2,-4,2,0,-2,-2,8,1],[4,2,0,0,0,1,-2,-4,-1,-3,-3,-1,-1,2,1,1,3,-5,-3,-1,-2,1,3,0,0,1,10]]
G=[ [x//2 for x in r[:26]] for r in G27[:26] ]
assert all(G27[i][j]==2*G[i][j] for i in range(26) for j in range(26))
V=[[-1,0,2,0,0,1,-2,-2,-2,0,1,0,1,0,1,-1,1,0,1,0,0,-1,-1,-2,0,-1],[1,-2,0,0,-2,1,0,0,0,0,1,0,1,0,-1,1,-1,-2,-1,0,0,-1,-1,0,0,-1],[1,0,0,0,2,-3,-4,0,-2,2,3,2,1,0,-1,1,-1,-2,-1,2,0,1,-1,0,0,-1],[3,-2,4,0,2,1,2,0,0,0,1,0,3,0,1,1,1,0,-1,0,2,-1,-1,0,2,-1]] # twice four nearest vectors in the original Bacher basis
def gv(v):return [sum(a*b for a,b in zip(r,v)) for r in G]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
assert [[dot(v,gv(w)) for w in V]for v in V]==[[4+8*(i==j)for j in range(4)]for i in range(4)]
assert all(all((i-j)%2==0 for i,j in zip(V[0],z))for z in V)
N=26
L=[[F(0)]*N for _ in G];D=[]
for j in range(N):
 D.append(F(G[j][j])-sum(D[k]*L[j][k]**2 for k in range(j)));assert D[-1]>0;L[j][j]=1
 for i in range(j+1,N):L[i][j]=(G[i][j]-sum(D[k]*L[j][k]*L[i][k]for k in range(j)))/D[j]
det=F(1)
for d in D:det*=d
assert det==3
# A proof-safe DFS: a partial sum of individually floored nonnegative LDL terms can only underestimate the actual norm.
scale=1<<20;coeff=[];den=[];top=[];bot=[]
for j in range(N):
 e=lcm(*(F(x).denominator for x in [L[i][j] for i in range(j,N)]))
 coeff.append([int(e*L[i][j]) for i in range(N)]);den.append(e);top.append(D[j].numerator*scale);bot.append(D[j].denominator*e*e)
def points(m,par=None):
 v=[0]*N
 def rec(k,left):
  if k<0:
   if dot(v,gv(v))==m:yield v.copy()
   return
  s=sum(coeff[k][i]*v[i] for i in range(k+1,N));e=den[k]
  B=isqrt(((left+1)*bot[k]-1)//top[k])
  lo=-((B+s)//e);hi=(B-s)//e
  step=1 if par is None else 2
  if par is not None:lo+=(par[k]-lo)%2
  for x in range(lo,hi+1,step):
   v[k]=x;yield from rec(k-1,left-(e*x+s)**2*top[k]//bot[k])
 yield from rec(N-1,m*scale)
def canon(p):
 q=[x%24 for x in p];q[0]+=24*(sum((x-y)//24 for x,y in zip(p,q))%2);return tuple(q)
def plus(p,q):return canon([x+y for x,y in zip(p,q)])
pr=[]
for i in range(N):
 a=[gv(v)[i] for v in V];pr.append([6*t-sum(a) for t in a]) # projected z-coordinate has denominator 24
zero=(0,0,0,0);qset={zero};queue=[zero]
for v in queue:
 for g in pr:
  t=plus(v,g)
  if t not in qset:qset.add(t);queue.append(t)
assert len(qset)==64
def qmin(c):
 # a representative within these limits has score <=3; any better one has |coord|<sqrt(3/2)
 return min(F(2*sum(x*x for x in p)+sum(p)**2,576)for n in product(range(-3,4),repeat=4) if sum(n)%2==0 for p in [[u+24*v for u,v in zip(c,n)]])
q={c:qmin(c)for c in queue}
pop=Counter()
for x in points(4):
 p=canon([sum(x[i]*pr[i][j]for i in range(N)) for j in range(4)]);pop[p]+=1
assert sum(pop.values())==117936
assert sorted(q[c] for c in queue if not pop[c])==[F(5,2)]*6+[3]*2
# Every provisional shell-6 upper bound is attained, so W is the exact required vertical metric.
for key,v in {"30,18,18,6":[1,-3,1,3,2,-1,-2,1,-1,-1,4,0,3,1,-2,1,1,-2,-2,0,0,-2,-1,-1,0,-2],"42,6,18,6":[1,-3,0,4,2,-2,-2,1,-1,-2,4,0,3,2,-2,0,1,-2,-1,0,0,-2,-1,-1,0,-2],"42,18,6,6":[1,-3,0,4,3,-2,-3,1,-2,-1,5,0,3,2,-2,0,2,-2,-1,0,0,-2,-1,-1,0,-2],"42,6,6,18":[2,-3,0,3,1,-2,-2,1,-1,-1,4,0,2,2,-2,1,1,-2,-1,1,0,-2,-1,-1,0,-2],"30,18,6,18":[-1,-2,0,3,1,-1,-4,0,-2,1,4,1,2,1,-2,1,0,-3,-2,0,-1,-1,-1,-1,0,-2],"36,12,12,12":[0,-2,0,2,0,-1,-3,0,-2,1,3,0,1,1,-2,1,0,-2,-1,0,-1,-1,-1,-1,0,-2],"30,6,18,18":[1,-3,-1,4,1,-3,-4,1,-1,1,4,1,2,1,-2,1,0,-3,-2,1,-1,-1,-1,-1,0,-2],"24,0,0,0":[1,-1,1,1,3,-1,-2,0,-2,-1,3,0,2,1,-1,0,1,-1,0,0,0,-1,-1,-1,0,-2]}.items():
 assert dot(v,gv(v))==6
 p=[sum(v[i]*pr[i][j] for i in range(26))for j in range(4)]
 qc=tuple(map(int,key.split(',')));assert canon(p)==qc and q[qc]==F(2*sum(t*t for t in p)+sum(p)**2,576)
gen=[[30,18,18,18],[18,6,18,6],[36,0,0,12],[36,12,0,0]]
def at(i):
 p=zero
 for a,g in zip([i%4,(i//4)%4,(i>>4)%2,i>>5],gen):
  for _ in range(a):p=plus(p,g)
 return p
assert len({at(i)for i in range(64)})==64
W=[max(0,q[c]-(0 if pop[c] else 2))for c in map(at,range(64))]
assert list(map(float,W))==[0,0.75,1,0.75,0.5,0.75,0.5,1.75,2,1.75,1,1.75,0.5,1.75,0.5,0.75,1,0.75,2,0.75,0.5,1.75,0.5,0.75,1,1.75,2,1.75,0.5,0.75,0.5,1.75,1,0.75,2,0.75,0.5,1.75,0.5,0.75,1,1.75,2,1.75,0.5,0.75,0.5,1.75,1,1.75,2,1.75,1.5,1.75,1.5,0.75,1,0.75,2,0.75,1.5,0.75,1.5,1.75]
# Direct D5 test: torsion basis b_i=e_i-e_5 (0<=i<4); b_4=2 e_5, each coefficient mod4.
add=lambda a,b:((a^b)^(((a&b)&341)<<1))&1023
def norm(i):
 u=[(i>>(2*j))%4 for j in range(5)]
 p=u[:4]+[2*u[4]-sum(u[:4])]
 b=[(x+2)//4 for x in p];r=[x-4*y for x,y in zip(p,b)]
 return F(sum(x*x for x in r)+(min(16-8*abs(x) for x in r)if sum(b)%2 else 0),8)
C=[norm(i)for i in range(1024)];halves=[i for i in range(1024)if not add(i,i)];rt={x:[i for i in range(1024)if add(i,i)==x]for x in halves}
bits=[0,2,8,10,16,18,24,26,32,34,40,42,48,50,56,58]
def good(h):return all(C[x]>=W[i]for x,i in zip(h,bits))
accepted=0; hits=0
for b,a in product(halves,repeat=2):
 v=[0,a,b,add(a,b)]
 if not good(v):continue
 for c in halves:
  s=v+[add(x,c)for x in v]
  if not good(s):continue
  for d in halves:
   h=s+[add(x,d)for x in s]
   if not good(h):continue
   accepted+=1
   va=[p for p in rt[a] if all(C[add(p,x)]>=W[i+1] for x,i in zip(h,bits))]
   vb=[p for p in rt[b] if all(C[add(p,x)]>=W[i+4] for x,i in zip(h,bits))]
   for p,t in product(va,vb):
    hits+=all(C[add(add(p,t),x)]>=W[i+5] for x,i in zip(h,bits))
assert accepted==480 and hits==0
print('27D full glue target: 64 labels; weighted costs',list(map(float,W)), 'half assignments',accepted,'full D5 assignments',hits)

