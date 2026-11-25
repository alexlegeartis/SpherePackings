"""Exact enumeration of Selling-conorm vertices (integer arithmetic)."""
from itertools import combinations
from fractions import Fraction as Q
from functools import lru_cache
E=list(combinations(range(4),2))
cuts=[[((m>>a)^(m>>b))&1 for a,b in E] for m in range(1,8)]
M=cuts+[[int(i==j) for j in range(6)] for i in range(6)]
trees=[]
for t in combinations(range(6),3):
 s={0}
 for _ in range(3):
  for i in t:
   a,b=E[i]
   if a in s or b in s:s.update((a,b))
 if len(s)==4:trees.append(t)
def determinant(A):
 A=list(map(list,A));n=len(A);sg=den=1
 if n==0:return 1
 for j in range(n-1):
  k=next((k for k in range(j,n) if A[k][j]),None)
  if k is None:return 0
  if k!=j:A[k],A[j]=A[j],A[k];sg=-sg
  for i in range(j+1,n):
   for h in range(j+1,n):
    x=A[j][j]*A[i][h]-A[i][j]*A[j][h]
    assert not x%den
    A[i][h]=x//den
  den=A[j][j]
 return sg*A[-1][-1]
mins=[Q(10**9)]*7
rhs0=[[8 if (i&f).bit_count()%2 else 12 for i in range(1,8)]+[0]*6 for f in range(1,8)]
for sel in combinations(range(13),6):
 mat=[M[s] for s in sel];d=determinant(mat)
 if not d:continue
 for j,rhs in enumerate(rhs0):
  b=[rhs[s] for s in sel]
  u=[determinant([[b[i] if h==k else mat[i][h] for h in range(6)] for i in range(6)]) for k in range(6)]
  if d<0:u=[-x for x in u]
  dd=abs(d)
  if any(sum(i*k for i,k in zip(r,u))<dd*y for r,y in zip(M,rhs)):continue
  det=Q(sum(u[i]*u[k]*u[l] for i,k,l in trees),dd**3);mins[j]=min(mins[j],det)
print('minimum squared volumes by functional:',[str(z) for z in mins])
assert mins==[432]*6+[400]
