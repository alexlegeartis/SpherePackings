"""Independent integer replay of the closed-box certificates and incompatible neighbor claims (Python standard library)."""
import json,itertools,sys
from pathlib import Path
from math import gcd,sqrt,pi,gamma
from collections import Counter
G=(0,1,3,2)
C=sorted({sum(G[int(p[(i+k)%5])*s%4]<<(2*i) for i in range(5)) for p in ('01112','03110','21310','21132') for k in range(5) for s in (1,3)})
def wt(x): return x.bit_count()
assert len(C)==40 and min(wt(x^y) for x,y in itertools.combinations(C,2))==4
assert len({wt(c)%2 for c in C})==1
B=set(a for a in range(1024) if min(wt(a^c) for c in C)==3)
ds=Counter(min(wt(a^c) for c in C) for a in range(1024));assert ds=={0:40,1:400,2:472,3:112}
assert len({wt(b)%2 for b in B})==1 and len({wt(c)%2 for c in C}|{wt(b)%2 for b in B})==2
def perm(x,P): return sum(1<<P[i] for i in range(10) if x>>i&1)
def automorphisms():
 cs=set(C)
 for s,k,t in itertools.product([1,-1],range(5),range(32)):
  p=[2*((k+s*i)%5)+(j^(t>>i&1)) for i in range(5) for j in range(2)]
  permc={perm(c,p) for c in C}
  for q in C:
   z=q^perm(C[0],p)
   if {u^z for u in permc}==cs:
    yield p,z
    break
AU=list(automorphisms());assert len(AU)==320
def orbit(a,i):
 return {(perm(a,p)^t,(p[i] if i>=0 else -1)) for p,t in AU}
def domain():
 return {(a,i) for a in range(1024) for i in (range(10) if a in B else [-1])}
D=domain()
def representatives():
 rem=D.copy();reps=[]
 while rem:
  a,i=min(rem);o=orbit(a,i); assert o<=rem;rem-=o;reps.append((a,i))
 return reps
reps=representatives()
if len(sys.argv)==2 and sys.argv[1]!='--check':
 large=json.load(open(sys.argv[1]));assert C==large['code']
 small={'scale':large['scale'],'trees':[[a,i,large['trees'][str(a)][max(0,i)]] for a,i in reps]}
 Path(__file__).with_name('certificate.json').write_text(json.dumps(small,separators=(',',':'))+'\n')
cp=json.loads(Path(__file__).with_name('certificate.json').read_text())
assert [(a,i) for a,i,T in cp['trees']]==reps
S=cp['scale'];assert S==4096
counts=Counter()
def replay(a,l,u,t):
 if isinstance(t,dict):
  j=t['j'];assert 0<=j<10 and (l[j]+u[j])%2==0 and l[j]<u[j]
  m=(l[j]+u[j])//2;s=u.copy();s[j]=m;r=l.copy();r[j]=m
  replay(a,l,s,t['left']);replay(a,r,u,t['right']);counts['split']+=1
 else:
  p={k:int(v) for k,v in t};assert p and len(p)==len(t) and min(p.values())>0 and min(p)>=0 and max(p)<40
  h=sum(p.values());val=0
  for j in range(10):
   d=sum(v for k,v in p.items() if (a^C[k])>>j&1)
   val+=max((h-d)*x*x+d*(S-x)**2 for x in (l[j],u[j]))
  assert 8*val<=21*h*S*S;counts['leaf']+=1
for a,i,T in cp['trees']:
 l=[0]*10;u=[S//2]*10
 if i>=0: l[i]=S//4
 replay(a,l,u,T)
# rank of the binary differences, checked without a linear algebra library:
orth=[x for x in range(1024) if all(wt(x&(c^C[0]))%2==0 for c in C)]
assert orth==[0,1023]
# Graph on coordinate pairs of nearest triples around each bad binary corner.
# Every vertex has at least six neighbors, and no independent triple exists.
for a in B:
 triples=[c^a for c in C if wt(c^a)==3]
 adj=[[any(t>>i&1 and t>>j&1 for t in triples) for j in range(10)] for i in range(10)]
 assert min(sum(row)-1 for row in adj)>=6
 assert all(adj[i][j] or adj[i][k] or adj[j][k] for i,j,k in itertools.combinations(range(10),3))
witness=[459,711,492,846,141,900,989,994,503,833,263,225,222,130,464,68,813,523,958,544,46,17,945,567,540,786,315,184,125,114,763,888]
assert len(set(witness))==32 and set(witness)<=B
F=C+[x+1024 for x in witness]
assert min(wt(a^b) for a,b in itertools.combinations(F,2))==4
assert len(F)==72
print('112 deep corners; 320 checked affine maps;',len(D),'boxes reduced to',len(reps),'roots;',dict(counts))
print('40+32 construction, center density 9/256, occupied',9/256*pi**5.5/gamma(6.5))
for k,delta,ref in ((1,1/2,9/256),(2,1/(2*sqrt(3)),1/27),(3,1/(4*sqrt(2)),9/256)):
 b=40/1024*(2/sqrt(11/8))**k*delta
 print(10+k, 'full-Best carrier upper',b,'incumbent',ref,'ratio',b/ref)

