"""Integer count for the sufficient-condition chart in NONINVARIANT.md.
The geometric input is the proved depth>=2 property of verify_fifth_complete.py.
"""
from fractions import Fraction

p=47
C={
    3:(p+1)*p**2,
    4:(p+1)*p**2*(p**2+p)+(p+1)*p,
    5:(p**6-1)//(p-1),
}
assert C=={3:106032,4:239210448,5:234330768}
def gauss(n,k):
    r=Fraction(1,1)
    for j in range(k):
        r*=Fraction(p**(n-j)-1,p**(k-j)-1)
    assert r.denominator==1
    return r.numerator
expected={
    5:5701015263809418421488,
    6:1167788319653328,
    7:106032,
}
for d in (5,4,3):
    exponent=10-d
    terms=[C[e]*gauss(4,d-e)*p**((6-e)*(d-e)) for e in range(3,d+1)]
    ans=sum(terms)
    assert ans==expected[exponent]<gauss(10,d)
    assert 24**46<16**46*p**exponent
    print('Index',f'{p}^{exponent}', 'exact size of short-witness chart',ans)
mixed6=(p+1)*p**3+(p+1)*p+1
mixed7_25_34=(p+1)*p**2+p+1
assert mixed6==4985761 and mixed7_25_34==106080
print('Mixed invariant index 47^6 all excluded:',mixed6)
print('Types (2,5) and (3,4) at invariant index 47^7 excluded:',mixed7_25_34)
