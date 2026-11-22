"""Independent coordinate/field/density cross-check. No external packages."""
from fractions import Fraction as F
from itertools import product
from math import factorial
from decimal import Decimal as D, localcontext
p=47
# (a,b) represents a + bω, with ω² + ω + 1 = 0.
def mul(u,v):
 a,b=u;c,d=v
 return ((a*c-b*d)%p, (a*d+b*c-b*d)%p)
def pow0(u,n):
 z=(1,0)
 for _ in range(n):z=mul(z,u)
 return z
assert all((u*u+u+1)%p for u in range(p))
t=(3,23);c=(4,8)
assert pow0(t,24)==c and mul(c,c)==(p-1,0)
assert sum(pow0(u,48)==(p-1,0) for u in product(range(p),repeat=2))==48
# Invert the specified simple-root coordinates by parity.
stock={}
for a,b,d,e in product(range(-4,5),repeat=4):
 if (a+b+d+e)%2 or not (0<a*a+b*b+d*d+e*e<18):continue
 v=(a,a+b,(a+b+d-e)//2,(a+b+d+e)//2)
 assert (v[0],-v[0]+v[1],-v[1]+v[2]+v[3],-v[2]+v[3])==(a,b,d,e)
 assert any(sum(x*y for x,y in zip(v,row))%47 for row in ((1,26,36,35),(0,33,35,37)))
 w=a*a+b*b+d*d+e*e
 stock[w]=1+stock.get(w,0)
assert sum(stock.values())==624
assert 33%47
cent=F(9,2)**48/(2*47**9)
old2=F(470782630154001057287684206740574981207676658334802588467200000000000000000,5104921439776057931517782239322133981907513)
assert (cent*10000)**2>10516**2*old2 and (cent*10000)**2<10517**2*old2
# Directed rational π enclosure (Machin).
def terms(t,n):return sum((F((-1)**j,(2*j+1)*t**(2*j+1))for j in range(n)),F(0))
lo=16*terms(5,64)-4*terms(239,65)
hi=16*terms(5,65)-4*terms(239,64)
lower=cent*lo**48/factorial(48);upper=cent*hi**48/factorial(48)
assert lower>F('5.9371e-22') and upper<F('6.25082581067615e-22')
with localcontext() as ctx:
 ctx.prec=48
 f=lambda x:D(x.numerator)/D(x.denominator)
 print('norm histogram',stock,'center',f(cent),'center squared',cent**2)
 print('occupancy',f(lower),f(upper),'relative',((f(cent**2/old2)).sqrt()-1)*100)
