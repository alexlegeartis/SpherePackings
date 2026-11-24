"""All 135 E8 isotropic parity labels have an orthogonal 8-frame."""
from itertools import combinations,product
from pathlib import Path

def run():
 def q(x): return (x&(x>>4)).bit_count()%2
 def B(x,y): return q(x^y)^q(x)^q(y)
 def dot(x,y):return sum(a*b for a,b in zip(x,y))
 data=list(map(lambda l:tuple(map(int,l.split())),Path(__file__).with_name('roots.in').read_text().splitlines()))
 assert len(data)==256
 roots=[x for x in range(256) if q(x)]
 assert len(roots)==120
 for x in roots:
  u=data[x]
  assert len(u)==8 and len({v%2 for v in u})==1 and sum(u)%4==0 and dot(u,u)==8
  for y in roots:
   z=dot(u,data[y])
   assert z%4==0 and z//4%2==B(x,y)
 # These roots span E8/2E8 because their known form is nondegenerate.
 S={0}
 for x in roots: S|={x^v for v in S}
 assert len(S)==256
 labels={tuple(B(x,y) for y in roots):x for x in range(256)}
 assert len(labels)==256
 R=[]
 for i in range(8):
  R += [tuple(4*a*(j==i) for j in range(8)) for a in (1,-1)]
 for J in combinations(range(8),4):
  for s in product((1,-1),repeat=4):
   v=[0]*8
   for i,a in zip(J,s):v[i]=2*a
   R.append(tuple(v))
 for i in range(8):
  for s in product((1,-1),repeat=8):
   v=tuple(a*(3 if j==i else 1) for j,a in enumerate(s))
   if sum(v)%4==0:R.append(v)
 assert len(set(R))==len(R)==2160
 F={}
 for v in R:
  assert len({x%2 for x in v})==1 and sum(v)%4==0 and dot(v,v)==16
  s=[dot(v,data[y]) for y in roots]
  assert all(z%4==0 for z in s)
  x=labels[tuple(z//4%2 for z in s)]
  assert x and not q(x)
  if tuple(-z for z in v)<v:continue
  F.setdefault(x,[]).append(v)
 assert len(F)==135
 for V in F.values():
  assert len(V)==8 and all(dot(a,b)==0 for a,b in combinations(V,2))
  # Correlation follows from Pythagoras in the basis of V:
  # given norm-2 root, integer pairings with the 8 frame vectors have
  # sum of squares 8; for a norm-4 vector, sum of squares 16.
  for R0,r0,z0 in (([data[y] for y in roots],8,4),(list(R)[:1],16,8)):
   for v in R0:
    av=[dot(v,y) for y in V]
    assert all(z%4==0 for z in av) and sum((z//4)**2 for z in av)==r0
    assert max(map(abs,av))>=z0
 print("PASS 120 anisotropic root words and 135 isotropic orthogonal 8-frames")
if __name__=='__main__':run()
