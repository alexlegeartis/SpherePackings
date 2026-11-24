from itertools import product
from fractions import Fraction as F
p=41
def mul(u,v):
    a,b=u;c,d=v
    return ((a*c-b*d)%p,(a*d+b*c-b*d)%p)
def powf(v,n):
    x=(1,0)
    for _ in range(n):x=mul(x,v)
    return x
def vertical():
    for x in product(range(-3,4),repeat=5):
        q=sum(a*a for a in x)
        if not q or q>14 or sum(x)%2:continue
        u=x[0];v=u+x[1];w=v+x[2]
        a=(w+x[3]-x[4])//2;b=a+x[4]
        assert (u,-u+v,-v+w,-w+a+b,-a+b)==x
        yield u,v,w,a,b
def ok(C,R):
    return all(any(sum(a*u for a,u in zip(z,c))%p for c in zip(*C)) for z in R)
def density_floor(mu,q):
    # Covering-independent exact centre radius power and rational pi lower.
    from math import factorial
    from decimal import Decimal,localcontext
    # Volume V_85 = 2^86*π^42*43!/86! (equivalently 2^43 π^42/85!!).
    pi=16*sum(F((-1)**i,(2*i+1)*5**(2*i+1))for i in range(80))-4*sum(F((-1)**i,(2*i+1)*239**(2*i+1))for i in range(81))
    sq=(F(2)**86*pi**42*F(factorial(43),factorial(86)))**2*(F(mu,4)**85)/(q*q)
    # V85 formula from Gamma(43.5): actually 2^86 π^42 * 43! / 86!.
    with localcontext() as c:
        c.prec=45
        x=(Decimal(sq.numerator)/Decimal(sq.denominator)).sqrt()
    assert F('3.5478756025421e-19')**2<sq
    print("certified occupied density >=",F('3.5478756025421e-19'),"actual approx",x)
if __name__=='__main__':
    # All nonzero twists and choices are literal in F41(ω); ω²+ω+1=0.
    W=(0,9)
    assert powf(W,42)==(p-1,0)
    c=powf(W,21)
    assert c==(9,0) and mul(c,c)==(p-1,0)
    R=tuple(vertical())
    assert len(R)==2330
    C=((1,0),(26,33),(36,35),(35,37),(17,6))
    assert C[0]==(1,0) and any(v for u,v in C)
    assert ok(C,R)
    print("t0",W,"rho(P)",c,"vertical color columns",C,"short vectors",len(R))
    density_floor(16,2*p**8)
