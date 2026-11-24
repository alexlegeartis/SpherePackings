from fractions import Fraction as F
from math import isqrt
from itertools import product
from anchored_replay import G as G10,q
def G12(p,c,s):
 m=(12+s)//2; n=(36+s)//2
 return ((10,0,p,c),(0,10,c,p+c),(p,c,12,m),(c,p+c,m,n))
def inverse(G):
 A=[[F(x)for x in t]+[F(i==j) for j in range(4)]for i,t in enumerate(G)]
 for i in range(4):
  if not A[i][i]:return None
  assert A[i][i]>0
  A[i]=[x/A[i][i]for x in A[i]]
  for j in range(4):
   if i!=j:
    h=A[j][i];A[j]=[x-h*y for x,y in zip(A[j],A[i])]
 return [t[4:] for t in A]
def coeffs(G,M):
 H=inverse(G)
 if H is None:return None
 I=[isqrt(int(M*H[i][i]))for i in range(4)]
 return product(*(range(-i,i+1) for i in I))
def s(G,t):
 u=(2*t[1]-t[0],2*t[0]+t[1],2*t[3]-t[2],2*t[2]+t[3])
 return sum(x*y*a for x,row in zip(t,G)for y,a in zip(u,row))
def admissible(type,p,c,sgn):
 G=(G10 if type==10 else G12 if type==12 else G14)(p,c,sgn)
 C=coeffs(G,16)
 if C is None:return type==10 and max(abs(p),abs(c))==10
 for t in C:
  if any(t):
   n=q(G,t)
   if n>16:continue
   b=abs(s(G,t))
   if n<10 or (n==10 and b not in (2,6,10)) or (n==12 and b not in (0,4,8,12,16)) or (n==14 and b not in (2,6,10,14,18,22)) or (n==16 and b not in (0,4,8,12,16,20,24,28)) or b*b>=5*n*n:
    return False
 return True
if __name__=='__main__':
 from joint import jointrows
 R,V,M,A,B=jointrows()
 for b in (A,B):
  j=10 if b is A else 12
  ix=[i for i,t in enumerate(b)if any(not admissible(j,*z)for z in t)]
  print('removed ',j,len(ix),ix)

def G14(p,c,s):
 m=(14+s)//2; n=(42+s)//2
 return ((10,0,p,c),(0,10,c,p+c),(p,c,14,m),(c,p+c,m,n))
