"""Exact rank-96 square-affine closure and determinant-normalized metric neighborhood."""
from fractions import Fraction as F
from pathlib import Path
import json

P=97
S=F(9848857802,10**9)
U=F(11209)
W=U**2
assert S*S>P and 5604**2-97*569**2==-1 and 5604+569*S<U
Q={i*i%P for i in range(1,P)}
J=[[0 if i==j else (1 if (i-j)%P in Q else -1) for j in range(P)] for i in range(P)]
assert all(sum(row)==0 for row in J)
assert all(sum(J[i][t]*J[t][j] for t in range(P))==P*(i==j)-1 for i in range(P) for j in range(P))
M=[[pow(i,k,P) for i in range(P)] for k in range(48)]

source=[]
for v in json.loads(Path(__file__).with_name('vectors.json').read_text())['vectors']:
    assert len(v)==P and all(type(x) is int for x in v) and sum(v)==0 and any(v)
    q=sum(x*x for x in v)
    g=sum(x*y*J[i][j] for i,x in enumerate(v) for j,y in enumerate(v))
    assert g<=0
    a,b=(q+F(g)/S)/2,(q-F(g)/S)/2
    assert a>0 and b>0
    m=next((k for k in range(1,48) if sum(u*y for u,y in zip(v,M[k]))%P),48)
    source.append((m,a,b))

def envelope(pairs):
    """Breakpoints for piecewise-linear lower envelope on w=u² >= 1."""
    E=[]
    for a,b in sorted(set(pairs),key=lambda e:(-e[0],e[1])):
        if E and E[-1][0]==a:
            continue  # same slope, worse intercept
        left=F(1)
        while E:
            ap,bp,lp=E[-1]
            left=(b-bp)/(ap-a)
            if left>lp:break
            E.pop()
        E.append((a,b,left if E else F(1)))
    return E

baseline=F(9,2)**96/(4*47**18)   # squared center density of MW92-D4
margin=F(543,1000)
# Jensen gap for generalized condition ratio 11/8.
x=F(3,8)
t=1-sum((-1)**i*x**i/(i+1) for i in range(10))
y=F(455,543)  # 998/543 - 1
assert 48*(t*t/2+t**3/3+t**4/4+t**5/(5*(1-t))) < sum((-1)**(i-1)*y**i/i for i in range(1,21))
for m in range(1,49):
    E=envelope([(a,b) for k,a,b in source if k>=m])
    prev=F(1); n=0
    for j,(a,b,lo) in enumerate(E):
        hi=min(W,E[j+1][2]) if j+1<len(E) else W
        if lo>=W:break
        assert lo==prev and lo<hi
        # On [sqrt(lo),sqrt(hi)] the convex upper-bound for the
        # chosen nonzero vector is a*u+b/u, so test endpoints.
        bound4=max((a*w+b)**2/w for w in (lo,hi))
        squared_density=(bound4/16)**48/P**(2*m-1)
        assert squared_density < margin**2 * baseline
        prev=hi; n+=1
    assert prev==W
    print('PASS L_%d: %d intervals, margin and arbitrary relative deformation bound verified'%(m,n))
