"""Exact rational jets, moment count, and density endpoint for the proof."""
from fractions import Fraction as F
from math import prod

N = 7


def series(*a, n=N):
    return list(map(F, a)) + [F(0)] * (n - len(a))


def add(x, y):
    return [a + b for a, b in zip(x, y)]


def scal(c, x):
    return [c * a for a in x]


def mul(x, y):
    return [sum(x[j] * y[i-j] for j in range(i+1)) for i in range(len(x))]


def power(x, k):
    n = len(x)
    if not x[0]:
        assert isinstance(k, int) and k >= 0
        out = series(1, n=n)
        for _ in range(k):
            out = mul(out, x)
        return out
    if x[0] != 1:
        assert isinstance(k, int) and k >= 0
        return scal(x[0]**k, power(scal(1/x[0],x), k))
    z, out, u = [F(0)] + x[1:], series(1, n=n), series(1, n=n)
    c = F(1)
    for i in range(1, n):
        c *= F(k-i+1) / i
        u = mul(u, z)
        out = add(out, scal(c, u))
    return out


def jet_linear(cols, a=(1, 0, 0)):
    k = len(cols)
    M = [[x[i] for x in cols] + [F(a[i])] for i in range(k)]
    for i in range(k):
        j = next(j for j in range(i, k) if M[j][i])
        M[i], M[j] = M[j], M[i]
        M[i] = scal(1/M[i][i], M[i])
        for j in range(k):
            if j != i:
                M[j] = add(M[j], scal(-M[j][i], M[i]))
    return tuple(row[-1] for row in M)


# Level 2, variable q=exp(pi i tau); indices 0,2,4 give the Q jet.
th3, th4 = series(1, 2, 0, 0, 2), series(1, -2, 0, 0, 2)
f = power(mul(th3, th4), 2)
u = add(power(th3, 4), power(th4, 4))
v = power(f, 2)

def hjet(x, y):
    return mul(power(f, x), power(u, y))[::2][:3]

for spec, expect, jet0 in [
    ([(10,2),(12,0)], -2400, ((4,-128,-384),(1,-96,4512))),
    ([(9,3),(11,1)], -864, ((8,0,-6912),(2,-128,3392))),
]:
    cols = [hjet(*p) for p in spec]
    assert tuple(map(tuple,cols)) == jet0
    X = jet_linear(cols)
    assert sum(c*a[2] for c,a in zip(X,cols)) == expect

for spec, cof, leadfactor, lead in [
    ([(8,4),(10,2),(12,0)], (F(75,2048),F(21,256),F(11,128)), 16**4, 2400),
    ([(7,5),(9,3),(11,1)], (F(27,4096),F(37,512),F(27,256)), 4*16**3, 108),
]:
    X = jet_linear([hjet(*p) for p in spec])
    assert X == cof and X[0]*leadfactor == lead

# Level 4, variable Q=exp(2 pi i tau). Under the scaled TST,
# A->C, B->-B, C->A, so the cusp coefficient is evaluated at A/B=-1.
A, B, C = series(0,16), series(1,-8,24), series(1,8,24)
for a,b,u0,j2 in [
    (9,2,4,-1632), (F(17,2),F(5,2),F(7,2),-1216),
    (F(35,4),F(9,4),F(15,4),-1416),
]:
    X = mul(power(B,a),power(C,b))
    cols = [mul(X,B), mul(X,A)]
    sol = jet_linear(cols)
    assert sol == (1,u0)
    assert sum(c*r[2] for c,r in zip(sol,cols)) == j2

for a,b,u0,v0,target,actual in [
    (9,1,5,F(83,8),102,2),
    (F(17,2),F(3,2),F(9,2),F(33,4),304,96),
    (F(33,4),F(7,4),F(17,4),F(233,32),516,108),
    (F(35,4),F(5,4),F(19,4),F(297,32),177,2),
]:
    X = mul(power(B,a),power(C,b))
    cols = [mul(X,mul(B,B)), mul(X,mul(A,B)), mul(X,mul(A,A))]
    sol = jet_linear(cols)
    assert sol == (1,u0,v0)
    assert all(sum(c*r[i] for c,r in zip(sol,cols)) == (i==0) for i in range(3))
    b4 = int(4*b)
    assert F(b4) == 4*b
    assert abs(1-u0+v0)*2**b4 == target and target > actual

# Ordinary weight-24 extremal theta and the degree-10 moment for v²=10.
def sigma(k,n):
    return sum(d**k for d in range(1,n+1) if n%d==0)
E4 = series(1, *(240*sigma(3,n) for n in range(1,N)))
E6 = series(1, *(-504*sigma(5,n) for n in range(1,N)))
D = scal(F(1,1728), add(power(E4,3),scal(-1,power(E6,2))))
cands = [power(E4,6),mul(power(E4,3),D),power(D,2)]
coeff = jet_linear(cands)
assert coeff == (1,-1440,125280)
num = sum(k*v[3] for k,v in zip(coeff,cands))
assert num == 52416000
p = series(1)
for i in range(5):
    p = mul(p,series(-i*i,1))
mom = [F(1)]
for j in range(1,6):
    mom.append(mom[-1]*F(60*(2*j-1),48+2*(j-1)))
count = num*sum(a*b for a,b in zip(p,mom))/prod(25-i*i for i in range(5))
assert count == 200 and count/2 == 100
assert F(11,16)**53*16**5 < 1
print("PASS: half and quarter jets, 11-design dot-5 moment, normal-density comparison")

