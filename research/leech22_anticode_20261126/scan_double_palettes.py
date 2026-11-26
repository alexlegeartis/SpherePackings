import os, sys, math, json, itertools, time
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
from functools import lru_cache
from itertools import combinations,product
from collections import Counter
T=time.monotonic()
p=lambda x:x.bit_count()
Gol=[]
for m in range(4096):
 h=0
 for j in range(12):
  if m>>j&1:h ^=0xae3<<j
 h|=(p(h)%2)<<23;Gol.append(h)
assert Counter(map(p,Gol))=={0:1,8:759,12:2576,16:759,24:1}
def shell4():
 X=np.zeros((196560,24),dtype=np.int8);k=0
 for i,j in combinations(range(24),2):
  for a,b in product([-4,4],repeat=2):
   X[k,i]=a;X[k,j]=b;k+=1
 signs=np.array([[2-4*(b>>i&1) for i in range(8)] for b in range(256) if p(b)%2==0],dtype=np.int8)
 for h in Gol:
  if p(h)==8:
   X[k:k+128,np.array([j for j in range(24) if h>>j&1])]=signs;k+=128
 for h in Gol:
  s=np.array([1-2*(h>>j&1) for j in range(24)],np.int8)
  X[k:k+24]=s
  X[k+np.arange(24),np.arange(24)] -=4*s;k+=24
 assert k==196560
 return X
X=shell4(); print('X',X.shape,flush=True)
def dprod(x,y):return int(np.dot(x.astype(np.int32),y.astype(np.int32)))//8
def mroot(i,j,s=1,t=-1):
 x=np.zeros(24,np.int16);x[i]=4*s;x[j]=4*t;return x
def odd(h=0,j=0,ty=5):
 x=np.array([1-2*(h>>i&1) for i in range(24)],np.int16);x[j]*=ty;return x
rng=np.random.default_rng(106)
norms=[]
for i in [0,1,2,3]:
 for j in range(i+1,6):
  for t in [-1,1]:norms.append(mroot(i,j,t,1))
Oct=[h for h in Gol if p(h)==8];Dod=[h for h in Gol if p(h)==12]
for j in range(5):
 for h in [0,*rng.choice(Oct,size=7,replace=False),*rng.choice(Dod,size=5,replace=False)]:
  for a in [5,-3]:norms.append(odd(int(h),j,a))
for j in range(14):
 for k in range(2):
  v=np.zeros(24,np.int16);h=Oct[j];J=[i for i in range(24) if h>>i&1]
  v[J]=2;v[J[0]]=-2;outside=[i for i in range(24) if i not in J]
  v[outside[k]]=4;norms.append(v)
for s in range(16):
 v=np.zeros(24,np.int16);ix=rng.choice(24,4,replace=False);v[ix]=4*rng.choice([-1,1],4);norms.append(v)
V=np.array(norms,dtype=np.int16)
root_ids=list(range(3))+[len(norms)-1]+[next(i for i,v in enumerate(norms) if dprod(v,v)==6)]
rng.shuffle(root_ids)
@lru_cache(maxsize=25)
def col(i):
 v=V[i]; out=X @ v
 assert np.all(out%8==0)
 return (out//8).astype(np.int16)
def atoms(a,b,c):
 D=a*c-b*b
 return sorted({((c*x-b*y)%D,(-b*x+a*y)%D) for x in range(D) for y in range(D)})
def geom(a,b,c):
 D=a*c-b*b
 labels=atoms(a,b,c); ids={h:j for j,h in enumerate(labels)};assert len(labels)==D
 smin=[10**9]*D
 R=max(a+abs(b),c+abs(b))
 for x in range(-R,R+1):
  for y in range(-R,R+1):
   i=ids[( (c*x-b*y)%D,(-b*x+a*y)%D)]
   smin[i]=min(smin[i],c*x*x-2*b*x*y+a*y*y)
 # parallelogram centered reps establish R sufficient? any smaller norm affects lower bounds only: using an upper bound on smin is always conservative.
 grp=[]
 for t in [0,1]:
  for k in labels:
   grp.append([ids[((h[0]-k[0])%D,(h[1]-k[1])%D)]+D*(t^s) for s in [0,1] for h in labels])
 return labels,ids,np.array(smin,np.int64),grp
geom=lru_cache(None)(geom)
@lru_cache(None)
def wanted(D,q):
 n=1+int(3*math.sqrt(D/12)*(44*D/(3*q))**11)
 def wins(k):return 12*k*k*(3*q)**22>9*D*(44*D)**22
 while n>1 and wins(n-1):n-=1
 while not wins(n):n+=1
 return n
@lru_cache(maxsize=4000)
def clique(adjs,k):
 # colored BB decision; translation already anchored at vertex 0
 n=len(adjs)
 def search(P,ls):
  m=k-len(ls)
  if m<=0:return ls
  if P.bit_count()<m:return None
  colors=[];order=[];R=P;cc=0
  while R:
   cc+=1;c=R
   while c:
    u=c&-c;j=u.bit_length()-1
    order.append(j);colors.append(cc);R^=u;c &= R & ~adjs[j]
  for j,cc in reversed(list(zip(order,colors))):
   if cc<m:return None
   ans=search(P&adjs[j],ls+(j,))
   if ans is not None:return ans
   P ^= 1<<j
  return None
 return search(adjs[0],(0,))
def test(D,F,G,use_double=True):
 N=D*(1+use_double)
 for q in sorted({min(x,16*D) for x in F[1:N] if x>0},reverse=True):
  k=wanted(D,q)
  if k>N:continue
  good=[x>=q for x in F[:N]];good[0]=False
  A=tuple(sum(1<<i for i,v in enumerate(row[:N]) if good[v]) for row in G[:N])
  if k>A[0].bit_count()+1:continue
  w=clique(A,k)
  if w is not None:return (q,k,w)
 return None
def cal(i,j,h):
 a,b,c,m=dprod(V[i],V[i]),dprod(V[i],V[j]),dprod(V[j],V[j]),dprod(V[h],V[h]);D=a*c-b*b
 u,v,w=col(i),col(j),col(h)
 if not all(np.any((u==p)&(v==q)&(w==r)) for p,q,r in [(1,0,0),(0,1,0),(0,0,1)]):return None
 labels,ids,smin,G=geom(a,b,c)
 k2=np.stack([((c*u-b*v)%D),((-b*u+a*v)%D)])
 lut=np.zeros((D,D),np.int16)
 for z,k in ids.items():lut[z]=k
 index=lut[k2[0],k2[1]]
 upper=np.full(2*D,16*D,dtype=np.int64)
 sq=(c*u*u-2*b*u*v+a*v*v).astype(np.int32)
 np.minimum.at(upper[:D],index,16*D-4*sq)
 np.minimum.at(upper[D:],index,16*D-4*sq+D*m-4*D*np.abs(w))
 upper[D]=min(upper[D],D*m)
 lower=np.minimum(upper, np.concatenate([24*D-4*smin,24*D-m*D-4*smin]))
 upper[0]=0;lower[0]=0
 # upper is capped at the same-fiber minimum 4: need norm4 witness perpendicular to pair
 assert np.any((u==0)&(v==0))
 key=(D,tuple(upper),tuple(lower),tuple(labels))
 return a,b,c,m,D,labels,G,list(map(int,lower)),list(map(int,upper)),key

def main():
 import hashlib
 seen={};results=[];attempt=0;P=[]
 for i in root_ids:
  for j in range(len(V)):
   if j==i:continue
   a,b,c=dprod(V[i],V[i]),dprod(V[i],V[j]),dprod(V[j],V[j])
   if abs(b)>min(a,c)//2:continue
   W=[z for z in range(len(V)) if dprod(V[z],V[i])==0 and dprod(V[z],V[j])==0]
   if not W:continue
   rng.shuffle(W)
   for m in [4,6,8]:
    WH=[z for z in W if dprod(V[z],V[z])==m]
    if WH:P.append((i,j,WH[0]))
 rng.shuffle(P)
 for i,j,h in P[:1200]:
  if time.monotonic()-T>480:break
  R=cal(i,j,h);attempt+=1
  if R is None:continue
  a,b,c,m,D,labels,G,lo,up,key=R
  if key in seen:continue
  seen[key]=1
  pos=test(D,lo,G); poss=test(D,up,G)
  rec={'gram':[a,b,c,m],'normals':[V[z].tolist() for z in (i,j,h)],'labels':labels,'low':lo,'high':up,'positive':pos,'undecided':poss}
  results.append(rec)
  if pos:print('WIN!',rec,flush=True);break
  if len(results)%15==0:print('models',len(results),'attempt',attempt,'possible',sum(r['undecided'] is not None for r in results),'elapsed',time.monotonic()-T,flush=True)
 print('FINAL',len(results),attempt,Counter((r['gram'][0],r['gram'][2],r['gram'][3]) for r in results),'unresolved',sum(r['undecided'] is not None for r in results),'positive',sum(r['positive'] is not None for r in results),'time',time.monotonic()-T,flush=True)
 out={'attempts':attempt,'geometry_count':len(results),'matrices':results}
 with open('rows.json','w') as f:json.dump(out,f,separators=(',',':'))
if __name__=='__main__':main()
