from fractions import Fraction as F
from math import prod
n=2**24
# theta_L = E4^3 -720 Delta; coefficients at q^2,q^3,q^4
def poly(a,b):
    out=[0]*5
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<5: out[i+j]+=x*y
    return out
E=[1]+[240*sum(d**3 for d in range(1,j+1) if j%d==0) for j in range(1,5)]
D=[0,1,0,0,0]
from math import comb
for j in range(1,5):
    D=poly(D,[(-1)**(i//j)*comb(24,i//j) if i%j==0 else 0 for i in range(5)])
s=[a-720*b for a,b in zip(poly(poly(E,E),E),D)]
assert s==[1,0,196560,16773120,398034000]
assert s[4]==48*(n-1-s[2]//2-s[3]//2)
rows=[(93150,94208,9200,0,2),(75900,97200,22356,1104,0),(65780,94208,32384,4096,92)]
for j,cs in enumerate(rows,2):
    assert sum(cs)==s[2]
    for h in range(1,6):
        assert sum(v*i**(2*h) for i,v in enumerate(cs))==F(s[2]*(8*j)**h*prod(range(1,2*h,2)),prod(range(24,24+2*h,2)))
k=s[2]//2
a,b,e=[k-sum(cs[1::2]) for cs in rows]
assert (k,a,b,e)==(98280,4072,-24,-24)
assert F(n*(-b),k-b)==4096
delta=F((k-a)*(k-b),n); assert delta==552
lam=a+b+delta; assert lam==4600
B=F(24**3*n*n,24**3+3*24*k+k*lam)
M=n*F(27,16)**12
assert B<M
assert F(16,27)**12<1 # singleton
assert F(8,9)**12<1 # no type-2 branch t>=1
# endpoint signs of the convex/concave inequalities used in the proof
for c in (F(0),F(1,2)):
    assert -24*c+40*c*c<=0
    assert -8-32*c+76*c*c<0
    assert F(1,2)+2*c-F(9,2)*c*c>0
    assert 6+2*c*c-2*c<=6
    assert 8+4*c*c-2*c<=8
print('PASS exact shell coefficients, design moments, spectrum and inequalities')
print('forbidden-support cap',B,'required motif count at hardest low endpoint',M)
