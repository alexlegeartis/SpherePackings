import sympy as s
from sympy import Rational as Q
for k in [2,3,4]:
 n=24*k
 N=2*k+5
 u=lambda i=1:[i]+[0]*(N-1)
 def add(a,b):return [i+j for i,j in zip(a,b)]
 def sm(c,a):return [c*x for x in a]
 def mul(a,b):return [sum(a[i]*b[j-i] for i in range(j+1)) for j in range(N)]
 def powr(a,b):
  p=u()
  for j in range(b):p=mul(p,a)
  return p
 def prod(*A):
  p=u()
  for a in A:p=mul(p,a)
  return p
 T=[int(any(i==j*(j+1)for j in range(10)))for i in range(N)]
 X=[0]+[16*x for x in powr(T,4)[:-1]]
 T3=[int(i==0)+2*int(i>0 and any(i==j*j for j in range(1,10)))for i in range(N)]
 T4=[c*(-1)**i for i,c in enumerate(T3)]
 Y=powr(T4,4);W=add(X,Y)
 e=2*k-1
 for j in [0,2,4]:
  d=2*k+1+j//2
  cols=[];sh=[]
  for i in range(d//2+1):
   a=prod(powr(mul(Y,W),e),powr(T3,2),powr(T4,2),powr(add(X,sm(2,Y)),d-2*i),powr(X,2*i))
   b=sm(4,prod(powr(mul(X,W),e),powr(T3,2),powr(T,2),powr(add(Y,sm(2,X)),d-2*i),powr(Y,2*i)))
   cols.append([a[2*l] for l in range(k+1)]);sh.append(b[e:e+2])
  M=s.Matrix(cols).T
  S=s.Matrix(sh).T
  if j==0:
   C=S*M.inv()*s.Matrix([1]+[0]*k)
   v0,v1=C
   assert tuple(C)=={2:(108,34928),3:(6272,2165248),4:(380632,136063200)}[k]
  else:
   C=S*M.nullspace()[0]
   assert C[0]
   if j==2:t2=C[1]/C[0]
   else:t4=C[1]/C[0]
   assert len(M.nullspace())==1 and C[1]/C[0]=={2:-20,4:4}[j]
 H=Q(4,3)*(2*k-Q(1,2)); d0=2*k-Q(1,2);d1=d0+1
 A2=-v0*H*d0/n*t2+v1*H*d1/n
 A4=v0*3*H**2*d0**2/((n+2)*(n+4))*t4 + 6*H*d1/(n+4)*A2-v1*3*H**2*d1**2/((n+2)*(n+4))
 c2=A2/H
 c4=A4/H**2
 C=(120*d1+24*d0)/(n+4)
 D=2*k
 r={2:Q(903,200),3:Q(3491,500),4:Q(9481,1000)}[k]
 assert c2>0 and c4>0 and C>0
 assert D-3*r/4>0 and (D-3*r/4)**2>r/12
 # The algebraic root lies above r, and the shallow branch lower.
 assert r>D and 3*r/4>D-1
 a,b=4*(r-d0),r
 assert 0<a<b<4*d0/3
 def f(h):
  t=min(d0+h/4-r,d0-3*h/4)
  assert t>=0
  return c4*h-C*v0*t*t-c2*(1+t)**2
 assert f(a)>0 and f(b)>0
 print(n,(v0,v1),'ratios',t2,t4,'positive endpoints',f(a),f(b))
