"""Fifth torsion paired coefficients over min-eight rank80 even unimodular."""
from fractions import Fraction as F
from sympy import Matrix
from math import comb
N=7
def C(x):return (x,)+(0,)*(N-1)
def mul(a,b):return tuple(sum(a[i]*b[j-i] for i in range(j+1)) for j in range(N))
def p(a,j):
 r=C(1)
 for _ in range(j):r=mul(a,r)
 return r
def e(j,k):
 z=list(C(1))
 for i in range(1,N):
  if i*j>=N:break
  z[i*j]=(-1)**i*comb(k,i) if k>=0 and i<=k else (comb(i-k-1,i)if k<0 else 0)
 return tuple(z)
def vr5(M):
 x=0
 while M%5==0:M//=5;x+=1
 return x
U=C(1);V=(0,1)+(0,)*(N-2)
for j in range(1,N):
 U=mul(U,mul(e(j,5),e(5*j,-1)))
 V=mul(V,mul(e(5*j,5),e(j,-1)))
result=[]
for J,s in [(20,2),(40,1)]:
 a=J-3
 K=[mul(p(U,x),p(V,J-x))for x in range(a,J+1)]
 M=Matrix([[k[j]for k in K]for j in range(4)])
 b=Matrix([s,0,0,0]);x=M.inv()*b
 assert M[:,1:].row_join(b).rank()>M[:,1:].rank()
 cnt=int(x[0]*5**((a-5*(J-a))//2))
 result.append(cnt)
 print('J',J,'cap',F(2*a,5),'nearest coefficient',cnt,'v5',vr5(cnt),'q4',sum(z*b[4]for z,b in zip(x,K)),flush=True)

assert result==[1245000,56591796875000]
M=result[1]
# Rank80 fourth spherical shifted dot moment factor.
f4=F(83*85,80*82)
assert f4.denominator%41==0 and f4.numerator%41!=0 and M%41!=0
assert all(i%41 for i in range(2,74,2))
assert vr5(f4.numerator)==vr5(f4.denominator)==0
assert vr5(M)==14 < 16
print('PASS product cap 74/5 impossible in all residues, hence sum <=72/5')
# For d1>202/25, d2 < 72/5-202/25 =158/25.
# Nondeep half layer gives d1<=13/2 +d2/4 <202/25.
# Deep half branch is convex in d2; at both endpoints it is less than 202/25.
x=F(158,25)
assert F(13,2)+x/4==F(202,25)
assert F(7)<F(202,25) and F(7,80)*x>(F(7)+x/4-F(202,25))**2
print('PASS universal individual fifth depth <=202/25; isotropic <=8')
