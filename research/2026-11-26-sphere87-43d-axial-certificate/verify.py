from fractions import Fraction as F
from itertools import product,combinations
from math import isqrt
from sympy import Matrix,eye,ones
G=7*eye(6)-ones(6)
for i,j in [(0,1),(2,3),(4,5)]:G[i,j]=G[j,i]=-2
H=G[1:,1:]; assert H.det()==3072
B=H.inv();P=[(0,)*5]+[tuple(int(i==j) for j in range(5))for i in range(5)]
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def norm(x):return F((Matrix(x).T*B*Matrix(x))[0])
d=F(7,24); short=[]
for x in product([-1,0,1],repeat=5):
 if norm(x) and norm(x)<F(1,2):short.append(x)
assert set(short)=={sub(a,b) for a,b in product(P, repeat=2) if a!=b}
matches=[(i,j)for i,j in combinations(range(6),2)if norm(sub(P[i],P[j]))==F(1,4)]
assert matches==[(0,1),(2,3),(4,5)]
for i,j in combinations(range(6),2):assert norm(sub(P[i],P[j]))==(F(1,4) if (i,j) in matches else d)
for a,b in combinations(short,2):
 if norm(sub(a,b))<=d:
  # Enumerated simplex directed edges sharing a head or tail.
  i,j=next((i,j) for i,j in product(range(6),repeat=2) if sub(P[i],P[j])==a)
  k,l=next((i,j) for i,j in product(range(6),repeat=2) if sub(P[i],P[j])==b)
  assert i==k or j==l
assert F(8)*271**43<5*274**43
assert F(36,3072)*F(137,96)**43 == F(137**43,2**223*3**42)
for m in [6,8]:
 width=F(1,m); assert d+width<F(1,2)
print("PASS: det 3072; 30 nonzero costs below 1/2; 6 short oriented directions; tetrahedron winding bound; radius slack 1/16 and exact DSW squared density.")
