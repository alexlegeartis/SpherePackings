# Python 3, no dependencies. Exhaustive disjoint-support certificate.
from functools import partial
N=18; M=(1<<N)-1
def qr(g):
    A=[0]
    for i in range(9):
        x=g<<i
        v=x|((x.bit_count()&1)<<17)
        A += [a^v for a in A]
    return sorted(A)
Q,R=qr(313),qr(471)
assert len(set(Q))==len(set(R))==512
assert set(Q)&set(R)=={0,M}
assert all(not ((c&d).bit_count()&1) for c in Q for d in R)
assert all(c.bit_count()>=6 for c in Q+R if c)
def four(E,forbid=()):
    F=tuple(x for x in E if x not in forbid)
    def go(i,u,n):
        if n==0:return ()
        for j in range(i,len(F)):
            v=F[j]
            if n*v.bit_count()>N-u.bit_count():break
            if u&v:continue
            k=go(j+1,u|v,n-1)
            if k is not None:return (v,)+k
        return None
    r=go(0,0,4)
    if r is not None:
        assert len(set(r))==4 and all(x not in forbid for x in r)
        assert all(not (a&b) for i,a in enumerate(r) for b in r[:i])
    return r
def bad(C):
    seen=bytearray(1<<N); out=set()
    for y in range(1<<N):
        if seen[y]:continue
        E=sorted((y^c for c in C),key=lambda x:(x.bit_count(),x))
        for x in E:seen[x]=1
        r=four(E);assert r is not None
        # All other antipodal pairs avoid this same valid quartet.
        for h in r:
            if four(E,(h,h^M)) is None:out.update((h,h^M))
    return out
badQ,badR=map(bad,(Q,R))
assert len(badQ)==len(badR)==344
C=sorted((set(Q)|set(R))-{0,M}); exc=set()
for y in badQ&badR:
    E=sorted((c^y for c in C),key=lambda x:(x.bit_count(),x))
    if four(E) is None:exc.add(y)
assert exc=={m for m in range(1<<N) if m.bit_count() in (0,1,17,18)}
print('all four-support certificates valid except the 38 analytically treated endpoints')
