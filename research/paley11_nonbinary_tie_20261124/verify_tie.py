from itertools import combinations
from fractions import Fraction as F
from pathlib import Path
R={i*i%11 for i in range(1,11)}
H=[[1]*12]+[[1]+[-1 if i==j else 1 if (i-j)%11 in R else -1 for j in range(11)]for i in range(11)]
assert all(sum(x*y for x,y in zip(a,b))==0 for a,b in combinations(H,2))
U=[sum((x==-1)<<i for i,x in enumerate(a))for a in H]
J={z for a,b in combinations(U,2)for z in [a^b,a^b^4095]}
B=sorted(x>>1 for x in J if x&1==0)
assert len(B)==66 and all(x.bit_count()==6 for x in B)
T=[(x,0)for x in B]
T.extend(tuple(map(int,x.split())) for x in Path(__file__).with_name('bank-best.txt').read_text().splitlines())
assert len(set(T))==72
X=[[(2*((v>>k)&1)+((s>>k)&1))%4 for k in range(11)]for v,s in T]
def norm(x,y):return sum(min((a-b)%4,(b-a)%4)**2 for a,b in zip(x,y))
assert min(norm(x,y)for x,y in combinations(X,2))==16
assert any(norm(x,y)%4==1 for x,y in combinations(X,2))
assert F(len(X),2**11)==F(9,256)
edges=[(i,j,u^v)for i,u in enumerate(B)for j,v in enumerate(B)if i<j and (u^v).bit_count()==4]
blocks={tuple(k for k in range(11)if w>>k&1)for i,j,w in edges}
assert len(edges)==990 and len(blocks)==330
assert len({s for b in blocks for s in combinations(b,2)})==55
for k in range(11):
 adj=[set()for x in B]
 for i,j,w in edges:
  if (w>>k)&1:adj[i].add(j);adj[j].add(i)
 C={0};todo=[0]
 for i in todo:
  for j in adj[i]-C:C.add(j);todo.append(j)
 assert len(C)==66
print("PASS distinct nonbinary 72-site 11D tie and infinitesimal strict rigidity inventory (66-site core)")
