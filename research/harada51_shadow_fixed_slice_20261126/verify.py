"""Exact verification of Harada51 geometry and all-radii 51/52/53 packing comparisons."""
from fractions import Fraction as F
from math import isqrt,sqrt,pi,gamma
import json
from pathlib import Path
def require(t):
 if not t: raise ArithmeticError("certificate did not verify")
R=[4,2,1,7,3,7,3,4,3,8,5,6,2,2,3,1,7,0,7,6,6,3,5,7,5]
D=25; n=51
def bal(x):return (x+4)%9-4
A=[[int(i==j) if j<D else 1 if j==D else bal(R[(j-D-1-i)%D]) for j in range(n)] for i in range(D)]
A.append([0]*D+[3]*(D+1))
A += [[9*int(i==j) for j in range(n)] for i in range(D+1,n)]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
v=[1]*25+[7]+[-1]*25     # integer numerator of the characteristic vector v, actual vector v/3
def isnum(w):
 k=w[25]-sum(w[:25])
 return len(w)==51 and k%3==0 and all((w[j+26]-sum(w[i]*R[(j-i)%25] for i in range(25))-k)%9==0 for j in range(25))
require(all(isnum(w) for w in A))
require(all(dot(x,y)%9==0 for x in A for y in A))
require(all(not A[i][j] for i in range(n) for j in range(i)))
require(A[0][0]==1 and A[25][25]==3 and A[50][50]==9)
require(isnum(v) and dot(v,v)==99 and dot(A[0],v)==9)
require(all((dot(x,x)-dot(x,v))%18==0 for x in A))
W=json.loads(Path(__file__).with_name('witness.json').read_text())
require(len(W)==6)
for m,w in enumerate(W):
 require(isnum(w) and dot(w,v)==9*m and dot(w,w)==9*[6,5,6,5,6,7][m])
H=[F(0)]+[F(m,11) for m in [54,62,46,50,52]]
# One-dimensional estimates: all same-color tails; then all different tails.
# Squared ratio to clean density at min6.
u=F(57,10)
small51=F(4,11)*121/u*(u/6)**51
hi51=F(4,11)/(u-H[2])*(u/6)**51
require(small51<1 and hi51<1 and F(4,11)/(6-H[2])==1)
# In 52D, low radii by the exact two-dimensional hexagonal bound used on each of 11 colors.
u=F(59,10)
small52=F(4,11)*121*F(4,3)/u**2*(u/6)**52
require(small52<F(49,50)**2)
# Higher radii. Lower bound pi and an upper left-end Riemann sum for each normalized disk cap.
N=10000
a=[0]*(N+1)
for j in range(N-1,-1,-1):
 t=N*N-j*j
 w=isqrt(t)
 if w*w<t:w+=1
 a[j]=a[j+1]+w # upper sum for integral sqrt(1-t^2)dt after division by N^2.
def low_area(l):
 R2=(l-H[5])/4
 def I(c):
  x2=(l-c)/(l-H[5]);require(0<x2<1)
  xlo=isqrt((N*N*x2.numerator)//x2.denominator)
  return F(a[xlo],N*N)
 return R2*(F(314159,100000)-4*(I(H[2])+I(H[1])))
worst=(F(0),None)
for j in range(5900,6000):
 l=F(j,1000);u=F(j+1,1000);B=low_area(l)
 require(B>0 and 4*(l-H[5])<l) # at most one near neighbor of each other color
 c=F(4,11)*(u/6)**52/B**2
 require(c<F(49,50)**2)
 if c>worst[0]:worst=(c,l)
# Exact endpoint oracles for the known 51–53D opening table.
for n0 in [51,52,53]:
 U=2 if n0%2 else 1
 if n0%2:
  for j in range(1,(n0-1)//2+1):U*=pi/(j+.5)
 else:
  for j in range(1,n0//2+1):U*=pi/j
 de=1.5**(n0/2)/2
 print(n0,"benchmark center:",de,"occupied:",U*de)
print("Six exact shortest projected-phase witnesses verified; minima numerators", [int(h*11) for h in H])
print("51D: small and high endpoint ratios",sqrt(small51),sqrt(hi51),1)
print("52D: cap-grid upper ratio",sqrt(worst[0]),"at interval starting",worst[1],"and simple low-range ratio",sqrt(small52))
print("53D: arbitrary tail upper ratio",sqrt(F(11,27)))


# Machin inequality used by both cap checks, with fully rational alternating Taylor tails.
def arct_bounds(n,last=6):
 s=sum(F((-1)**k,(2*k+1)*n**(2*k+1)) for k in range(last+1))
 return s-F(1,(2*last+3)*n**(2*last+3)),s
m1,m2=arct_bounds(5);n1,n2=arct_bounds(239)
require(16*m1-4*n2>F(314159,100000))
