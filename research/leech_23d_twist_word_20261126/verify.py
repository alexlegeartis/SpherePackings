"""Exact verification from this sibling checkout; no solver or floating point."""
from pathlib import Path
import runpy
P = Path(__file__).resolve().parent
old = runpy.run_path(str(P.parent/'leech_23d_projection_optimality_20261125'/'gluing72.py'))
D, dot, member, v, w = (old[k] for k in ('D','dot','member','v','w'))
seen = set()
for line in (P/'vectors.txt').read_text().splitlines():
    a,b,*x = map(int,line.split())
    assert len(x)==24 and member(x)
    A,rem=divmod(dot(x,v),8); assert not rem and A%6==a
    B,rem=divmod(dot(x,w),8); assert not rem and B%12==b
    assert 3*dot(x,x)//2-2*A*A-B*B==48-D(a,b)
    assert (a,b)!=(0,0) and (a,b) not in seen
    seen.add((a,b))
assert len(seen)==71 and old['H'][0,0]>0
G=[(a,b) for a in range(6) for b in range(12)]
def cap(h):
    adj=[sum(1<<j for j,(u,v) in enumerate(G) if D((u-a)%6,(v-b)%12)<h) for a,b in G]
    best=0
    def expand(P,n):
        nonlocal best
        best=max(best,n)
        while n+P.bit_count()>best:
            x=P&-P;P-=x;expand(P&adj[x.bit_length()-1],n+1)
    expand(adj[0],1)
    return best
assert [cap(h) for h in [3,6,12]]==[2,4,8]
for h,m,ends in [(3,2,[0,2]),(6,4,[2,4]),(12,8,[4,9]),(48,72,[9])]:
    for b in ends:
        gap=4*(h-b)*46**23-m*m*(48-b)**23
        assert gap==0 if b==2 else gap>0
print('71 exact distances, three clique numbers, all continuum endpoint bounds verified.')
