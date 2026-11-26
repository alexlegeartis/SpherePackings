"""Standalone exact separation and density certificate; standard library only."""
from fractions import Fraction
from itertools import combinations,product
rows=[0]+[sum(1<<(((i+j)%11)+1) for i in (0,1,3,4,5,9)) for j in range(11)]
C={w for a,b in combinations(rows,2) for w in (a^b,a^b^4095)}
assert len(C)==132 and all(c.bit_count()==6 for c in C)
C.update(w for i in range(0,12,2) for w in (3<<i,(3<<i)^4095))
assert len(C)==144 and min((a^b).bit_count() for a,b in combinations(C,2))==4
cosets=[[] for _ in range(32)]
for x in range(4096):
    syn=0
    for i in range(12):
        if (x>>i)&1: syn^=16+i
    cosets[syn].append(x)
assert all(len(a)==128 for a in cosets)
alpha,beta,gamma=28,29,30
def label(m,n):
    c=(m-n)%3
    if not c: return 1,(0,1,2)[n%3]
    return 0,((m%2)*(beta^gamma))^((n%2)*(alpha^gamma))^((c==2)*(alpha^beta^gamma))
def words(t,p):
    return [a^p for a in C] if t else cosets[p]
seen={}
for m,n in product(range(6),repeat=2):
    assert label(m,n)==label(m+6,n)==label(m,n+6)
    for dm,dn in product(range(-2,3),repeat=2):
        z=dm*dm+dm*dn+dn*dn
        if z>=4: continue
        t,p=label(m,n);s,q=label(m+dm,n+dn)
        key=(z,t,p,s,q)
        if key in seen: continue
        if s!=t: horizontal=3
        else:
            horizontal=min((x^y).bit_count() for x in words(t,p) for y in words(s,q)
                           if z!=0 or x!=y)
        seen[key]=z+horizontal
# Shifts outside this offset box cannot be shorter: Q≥(dm²+dn²)/2.
assert min(seen.values())==4 and len(seen)==113
count=sum(len(words(*label(m,n))) for m,n in product(range(6),repeat=2))
assert count==4800
coefficient=Fraction(count,2**12*18)
assert coefficient==Fraction(25,384) and coefficient/Fraction(1,16)==Fraction(25,24)
print("PASS: 4800 sites; exact min squared 4; center density 25/(384 sqrt(3)); ratio 25/24.")
