"""Verify the two 72-color gluings by exact signed minimal-shell arithmetic."""
from collections import Counter

v = [5] + [1] * 23
w = [0,-2,0,-2,-2,4,4,0,4,-2,0,2,0,0,0,0,-4,0,2,0,-2,0,0,-2]
G = sum(1 << j for j in (11,9,7,6,5,1,0))
C = set()
for k in range(4096):
    c = 0
    for j in range(12):
        if k >> j & 1:
            c ^= G << j
    C.add(c | ((c.bit_count() % 2) << 23))
assert len(C) == 4096
def dot(x, y):
    return sum(a*b for a,b in zip(x,y))
def member(x):
    p = x[0] % 2
    return (all((a-p)%2 == 0 for a in x)
            and (sum(x)-4*p)%8 == 0
            and sum(1 << i for i,a in enumerate(x) if (a-p)%4 == 2) in C)
assert member(v) and member(w) and (dot(v,v),dot(v,w),dot(w,w)) == (48,0,96)

# Three coordinate types of Leech norm-four vectors, in numerator units.
H = Counter()
def add(x):
    a,b = dot(v,x),dot(w,x)
    assert dot(x,x) == 32 and a%8 == b%8 == 0
    H[a//8,b//8] += 1
for i in range(24):
    for j in range(i):
        for a in (-4,4):
            for b in (-4,4):
                x = [0]*24
                x[i],x[j] = a,b
                add(x)
for c in C:
    I = [i for i in range(24) if c >> i & 1]
    if len(I) == 8:
        for k in range(256):
            if k.bit_count()%2 == 0:
                x = [0]*24
                for j,i in enumerate(I):
                    x[i] = -2 if k >> j & 1 else 2
                add(x)
    for j in range(24):
        x = [-1 if c >> i & 1 else 1 for i in range(24)]
        x[j] *= -3
        add(x)
assert sum(H.values()) == 196560 and H[1,0] > 0 and H[0,1] > 0

def D(i,j):
    if i == j == 0:
        return 48   # equal horizontal coset; excludes identical sites in the proof
    A = [a for a in range(-3, 4) if a % 6 == i]
    B = [b for b in range(-6, 7) if b % 12 == j]
    return 2*A[0]**2 + B[0]**2 - (0 if any(H[a,b] for a in A for b in B) else 24)

assert D(0,5) == 1
assert all(D(i,j) == D(i,5*j%12) for j in range(12) for i in range(6))
# Intralayer differences within one coset cost at least 4. Between layer
# numbers 0..6, check all difference colors; beyond six height alone suffices.
for advance in (1,5):
    assert min(48 - D(i % 6,advance*k % 12) + k*k
               for k in range(7) for i in (-1,0,1) if k or i) == 46
print("196560 shell sites, two primitive axes, exact minimum 46/12 in both packings")
