import sympy as s
from fractions import Fraction
from math import floor
q=s.Symbol('q')
def trunc(f):return s.series(f,q,0,4).removeO().expand()
def mul(a,b):return trunc(a*b)
A=1+6*q+6*q**3
U=s.Integer(1)
C=q
for n in range(1,4):
 U=mul(U,trunc((1-q**n)**3/(1-q**(3*n))))
 C=mul(C,trunc((1-q**(3*n))**9/(1-q**n)**3))
B=trunc(U**3)
print('A=',A,'U=',U,'B=',B,'C=',C)
def run(n,r):
 k=n//2-r;a=k%3; m=k//3; nrows=(n//24+1)
 D={}
 for i in range(m+1):
  j=m-i
  z=trunc(U**r*A**a*B**i*C**j)
  D[i]=[z.coeff(q,h)for h in range(nrows)]
 v=s.Matrix([1]+[0]*(nrows-1))
 for d in range(m,-1,-1):
  I=list(range(d,m+1))
  T=s.Matrix([D[i]for i in I]).T
  if T.row_join(v).rank()>T.rank(): continue
  sol=T.gauss_jordan_solve(v)[0]
  i=d;j=m-d
  e=(r-a+3*i-9*j)//2
  assert (r-a+3*i-9*j)%2==0
  count=sol[0]*s.Rational(3)**e
  print(n,r,'max possible d',d,'distance',Fraction(2*d,3)+Fraction(2*r,9),'mat',T.tolist(),'coeff',list(sol),'lead',count)
  return d,count
assert __debug__
expected={(48,0):(6,1836),(48,1):(5,105),(48,2):(5,450),
          (72,0):(9,105084),(72,1):(8,6060),(72,2):(8,25632)}
for (n,r),expect in expected.items():
 assert run(n,r)==expect
print("PASS: all six bounds and conditional equality-shell coefficients")
