"""Replay all rational inequalities used in the trace-twist proof; Python 3."""
from fractions import Fraction as F

def poly(roots):
    p = [1]
    for r in roots:
        a = [0] * (len(p) + 1)
        for j, x in enumerate(p):
            a[j] -= r*x
            a[j+1] += x
        p = a
    return p

def moments(D,n,m):
    return [1,0,F(m*m*D,n//2+1),0,F(3*m**4*D*D,(n//2+1)*(n//2+3)),0]

def reach(D,m,h):
    return F(2*m*h,h*h+D*m*m)

# Left endpoint then two adjacent grid pairs; the polynomial is nonnegative on the grid.
for D,n,m,h,step,r,E in [
    (2,48,6,3,1,[-2,0,1,2,3],F(-1152,25)),
    (2,72,8,4,1,[-3,0,1,3,4],F(-7040,481)),
    (5,72,8,6,2,[-4,0,2,4,6],F(-286720,481)),
    (5,96,10,6,2,[-4,0,2,4,6],F(-912000,833)),
]:
    assert sorted(r) == r and r[0] == -h+step
    assert r[2]-r[1] == r[4]-r[3] == step
    assert all(x % step == 0 for x in r)
    assert sum(p*u for p,u in zip(poly(r),moments(D,n,m))) == E < 0

for D,n,m in [(5,24,4),(5,48,6),(3,48,6)]:
    assert moments(D,n,m)[2] > 4  # cubic threshold on even grid
assert moments(2,48,6)[4] > 16  # even-grid refinement

for m,h in [(4,4),(6,4),(8,6),(10,6)]:
    assert m-F(h,5) > 0 and F(1,5)<reach(5,m,h)
    assert ((m-2)**2)*F(5,4) < m*m
assert 4*4*2 < 6*6

for n,m,h,b in [(48,6,3,8),(72,8,4,10)]:
    assert F(4,9) == reach(2,m,h)
    assert F(6*b-m,9) <= F(7*m,9)
    assert (6*b-m==7*m) == (n==48)
    assert F(b,2)**2 < F(m*m,2)
    assert 1-2*F(4,9)**2 == F(7,9)**2
assert F(4,9)<reach(2,6,4)

assert F(1,3)<reach(3,6,4)
assert F(14,3)**2<36*F(2,3)
assert F(4,2)<F(6,2)

# 96D rank-minimum 8 branch.
assert F(19,5)**2*5 < F(8619,1000)**2
assert F(8,1)**2/F(4,5) == 4**2*5
assert (5*10+3*(-14))//2 == 4 and (3*10-14)//2 == 8
assert (3*10-14)%4 == 0 and (3*10-10)%4 == 0
assert 10+F(-10,5) == 8
assert 12+F(16-3*12,5) == 8
print("PASS: moment witnesses, full-cone joins, and conditional 96D endpoint")
