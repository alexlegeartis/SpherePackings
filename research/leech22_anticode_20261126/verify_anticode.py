from collections import Counter
raw=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[2,-2,0,0,0,-2,2,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,2],[0,2,-2,0,0,0,2,2,-2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,2],[2,0,-2,0,0,0,2,0,-2,2,0,0,0,0,2,0,2,0,0,0,0,0,0,2],[-4,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[4,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,-1,1,-1,1,-1,3,1,-1,1,1,1,-1,1,-1,-1,1,-1,-1,1,-1,1,1,1],[1,1,-1,1,-1,-1,3,1,-1,1,1,1,-1,1,-1,1,1,1,-1,1,1,1,1,1],[1,-1,-1,-1,-1,-1,3,1,-1,-1,1,1,1,1,1,-1,1,1,-1,-1,-1,1,1,1],[-1,1,1,-1,-1,1,3,1,-1,1,1,1,1,1,-1,-1,1,1,1,-1,1,1,1,1],[1,-1,1,-1,-1,1,3,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,1,1,1],[-1,1,-1,1,1,1,3,1,-1,1,-1,1,1,1,1,-1,1,1,-1,1,-1,1,1,1]]
G=set()
for m in range(4096):
    c=0
    for j in range(12):
        if (m>>j)&1: c ^= 0xae3<<j
    G.add(c|((c.bit_count()%2)<<23))
assert Counter(c.bit_count() for c in G)=={0:1,8:759,12:2576,16:759,24:1}
w=[]
for x in raw:
    m=x[0]%2
    assert len(x)==24 and all((u-m)%2==0 for u in x)
    assert sum((((u-m)//2)%2)<<i for i,u in enumerate(x)) in G
    assert (sum(x)-4*m)%8==0
    w.append((sum(x[:3]),*x[3:]))
assert len(set(w))==12
# In these coordinates the physical squared distance is d24/24.
def d24(x,y):
    return (x[0]-y[0])**2 + 3*sum((a-b)**2 for a,b in zip(x[1:],y[1:]))
def profile(n):
    return Counter(d24(w[i],w[j]) for i in range(n) for j in range(i))
assert profile(4)=={72:6}
assert max(profile(12))==88 and sum(profile(12).values())==66
assert 4*9**11 < 11**11 and 12**11 < 3*11**11
print(profile(4),profile(12))
