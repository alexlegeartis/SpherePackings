"""Exact syndrome distances in E7^5/L35; independent 37D lift check."""
from itertools import product
from collections import Counter
# A different degree-six polynomial from the primary basis.
MOD=0b1100001
def mul(x,y):
 r=0
 while y:
  if y&1:r^=x
  x<<=1
  if x&64:x^=MOD
  y>>=1
 return r
def pw(x,n):
 p=1
 while n:
  if n&1:p=mul(p,x)
  x=mul(x,x);n>>=1
 return p
assert all(pw(x,63)==1 for x in range(1,64))
F8=[i for i in range(64)if pw(i,8)==i];assert len(F8)==8
def q(x):
 z=pw(x,9)
 return z^pw(z,2)^pw(z,4)
assert {q(x)for x in range(64)}=={0,1}
def cost(x,e):
 return 6*e if x==0 else 4-2*(q(x)^e)
D={0:0}
for a in F8[:5]:
 X={}
 options=[(x|(mul(a,x)<<6)|(e<<12),cost(x,e))for x,e in product(range(64),range(2))]
 for b,x in D.items():
  for s,y in options:
   ell=b^s
   if x+y<X.get(ell,1<<10):X[ell]=x+y
 D=X
assert len(D)==8192 and max(D.values())==6
for y in range(1,64):
 assert D[(y<<6)|(1<<12)]==6
 assert D[y<<6]==4
print('PASS all 8192 E7^5/L35 phase classes',sorted(Counter(D.values()).items()))
print('PASS 37D order-two lift: phase depth exactly 6; shifted by doubled order-four 36D phase depth exactly 4')
