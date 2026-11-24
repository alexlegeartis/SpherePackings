from fractions import Fraction
import sympy as s
for k in range(1,5):
 n=2*k+3
 I=lambda v:[v]+[0]*(n-1)
 def mul(a,b):return [sum(a[i]*b[j-i] for i in range(j+1)) for j in range(n)]
 def power(a,b):
  c=I(1)
  for _ in range(b):c=mul(c,a)
  return c
 x=power([int(any(i==j*(j+1) for j in range(10)))for i in range(n)],4)
 X=[0]+[16*t for t in x[:-1]]
 Y4=I(1);Y3=I(1)
 for a in range(1,n):
  if int(a**.5)**2==a:Y3[a]=2;Y4[a]=2*(-1)**int(a**.5)
 Y=power(Y4,4);Z=[i+2*j for i,j in zip(X,Y)];XY=[i+j for i,j in zip(X,Y)]
 pre=power(mul(Y,XY),2*k)
 def mat(d,extra):
  F=[mul(mul(pre,extra),mul(power(Z,d-2*j),power(X,2*j)))for j in range(d//2+1)]
  return s.Matrix([[f[2*j] for f in F]for j in range(k+1)])
 b=s.Matrix([1]+[0]*k)
 E=mat(2*k,I(1));E2=mat(2*k+1,I(1));O=mat(2*k-1,mul(power(Y3,2),power(Y4,2)))
 x=E.inv()*b
 print('rank',24*k,'iso minimum upper',2*k,'nearest count',(16**(2*k)*sum(x)), 'harmonic invert', E2.det()!=0, 'odd impossible',O.row_join(b).rank()>O.rank(), 'signed first shell', sum(a*c for a,c in zip(x,[mul(pre,mul(power(Z,2*k-2*j),power(X,2*j)))[2*k+2]for j in range(k+1)])))
