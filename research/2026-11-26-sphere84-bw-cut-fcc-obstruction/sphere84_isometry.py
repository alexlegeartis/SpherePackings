"""Normalize binary direction masks; check all integral Z^3 seed triples, including lifted-code isometry."""
from itertools import combinations
pc=int.bit_count
cols=[[8+(i%8) for i in range(16)], list(range(16)), [1,2,3,3,5,6,7,7,9,10,11,11,13,14,15,15]]
oldg=(254,855,1389)
def col19(c,g):return [v+sum(((z>>i)&1)<<(4+j) for j,z in enumerate(g)) for i,v in enumerate(c)]+[16,32,64]
base=col19(cols[1],oldg)
def iso(target):
 f=[base.count(i) for i in range(128)]; h=[target.count(i) for i in range(128)]
 def go(s,t):
  if len(s)==128:return True
  z=max((x for x in base if x not in s), key=lambda x:sum(f[x^y] for y in s))
  for w in target:
   if w not in t and all(f[z^x]==h[w^y] for x,y in zip(s,t)) and go(s+[z^x for x in s],t+[w^y for y in t]):return True
  return False
 return go([0],[0])
for k,c in enumerate(cols):
 U=[sum((pc(a&v)%2)<<j for j,v in enumerate(c)) for a in range(16)]
 D=[];deep=[]
 for m in range(65536):
  t=[m^u for u in U];d=min(map(pc,t));D.append(d)
  if min(t)==m and d>=7 and pc(m)%2:deep.append(m)
 T=[(x,y,z) for x,y,z in combinations(deep,3) if min(D[x^y],D[x^z],D[y^z])>=6 and D[x^y^z]>=5]
 assert len(T)==[0,448,448][k]
 for g in T:assert iso(col19(c,g))
 if k==2: assert (1367,1643,16761) in T
 print(k, sorted(map(pc,U)), len(deep),len(T))
print('All 896 seed triples give isometric 19D integral lifts.')

