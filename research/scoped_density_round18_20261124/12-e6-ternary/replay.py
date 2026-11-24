from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import isqrt

n = 729
V = list(product(range(3), repeat=6))
index = {x:i for i,x in enumerate(V)}
E = ((0,1),(1,2),(2,3),(3,4),(2,5))
G = [[2 if i==j else -int((i,j) in E or (j,i) in E) for j in range(6)] for i in range(6)]
def Q(x):
    return 2*(sum(t*t for t in x)-sum(x[i]*x[j] for i,j in E))
# Verify squared-norm coordinate bounds by supplying the exact 3*Gram inverse.
Ginv3 = [
    [4,5,6,4,2,3], [5,10,12,8,4,6], [6,12,18,12,6,9],
    [4,8,12,10,5,6], [2,4,6,5,4,3], [3,6,9,6,3,6]
]
assert all(sum(G[i][k]*Ginv3[k][j] for k in range(6)) == 3*(i==j)
           for i in range(6) for j in range(6))
B=[isqrt(12*Ginv3[i][i]//3) for i in range(6)]
r=(2,1,0,2,1,0)
shift=lambda x,k:tuple((a+k*b)%3 for a,b in zip(x,r))
d=[100]*n
for x in product(*(range(-a,a+1) for a in B)):
    norm=Q(x)
    if norm<=12:
        i=index[tuple(t%3 for t in x)]
        d[i]=min(d[i],norm)
rad={x for x in V if all(sum(G[i][j]*x[j] for j in range(6))%3==0
                         for i in range(6))}
assert rad == {shift(V[0],k) for k in range(3)}
assert {V[i] for i,v in enumerate(d) if v>8}==rad-{V[0]}
assert max(d)==12 and {d[index[u]] for u in rad-{V[0]}}=={12}
for x in V:
    if x not in rad:
        assert min(d[index[shift(x,k)]] for k in range(3)) == {0:6,1:2,2:4}[Q(x)//2 %3]

def dif(i,j):
    return index[tuple((a-b)%3 for a,b in zip(V[i],V[j]))]

S = (0,1,28,37,40,428)
U0 = (104,109,118,121,145,146,148,149,157,158)
U1 = (162,163,216,234,240,648,649,702,720,726)
U2 = (347,352,361,364,388,389,391,392,400,401)
U = U0 + U1 + U2
assert len(set(S+U))==36
assert all(d[dif(i,j)]==2 for i,j in combinations(S,2))
assert all(sum(d[dif(i,j)]==2 for j in S)>=4 for i in U)
assert all(d[dif(i,j)] in (2,4) for K in (U0,U1,U2)
           for i,j in combinations(S[1:]+K,2))

# Point orbits under the six exact simple reflections; no external orbit data.
def reflect(x,j):
    y=list(x)
    y[j]=(y[j]-sum(x[i]*G[i][j] for i in range(6)))%3
    return tuple(y)
R = [[index[reflect(x,j)] for x in V] for j in range(6)]
orbits, oid = [], [-1]*n
for i in range(n):
    if oid[i]>=0:
        continue
    S0={i}; todo=[i]
    while todo:
        x=todo.pop()
        for p in R:
            y=p[x]
            if y not in S0:
                S0.add(y); todo.append(y)
    for j in S0:
        oid[j]=len(orbits)
    orbits.append(sorted(S0))
sizes=list(map(len,orbits))
assert sizes == [1,72,270,240,72,72,1,1]
# H is twice the real character sum. Real part suffices even on split
# Weyl orbits: Fourier transforms of symmetrized autocorrelations are real.
H=sorted({tuple(3*sum(sum(a*b for a,b in zip(ch,V[i]))%3==0 for i in O)-len(O)
                for O in orbits)
          for ch in V})
assert len(H)==6

def stats(A,B):
    h=Counter(oid[dif(i,j)] for i in A for j in B)
    return [h[j] for j in range(len(orbits))]
ss,uu,su = stats(S,S),stats(U,U),stats(S,U)

# Signs of these exact dual multipliers are checked, not presumed.
a = [F(-118,119), F(2309,24276), F(253,1734), F(-361,48552),
     F(-1369,48552), F(1410,2023), F(-1535,163863), F(636,2023)]
b = [F(0),F(-2735,24276),F(-6935,163863),F(-24505,1310904),F(0),F(0)]
assert all(w<=0 for w in b)
for x,y in product(range(7),range(31)):
    polynomial = sum(c*t for c,t in zip(a,(x,y,x*x,y*y,x*y,1)))
    assert polynomial <= int(x==0 and y<=17)
for j in range(len(orbits)):
    lhs = a[2]*ss[j]+a[3]*uu[j]+a[4]*su[j]
    rhs = a[6]*sizes[j]+a[7]*(j==0)-sum(w*h[j] for h,w in zip(H,b))
    assert lhs>=rhs
q=F(127,n)
lo=sum(c*t for c,t in zip(a,(6*q,30*q,0,0,0,1,n*q*q,q)))
assert lo>0
assert 3*n-2*128==1931
volume=3*F(9,5)**3*F(18,5)**3
assert 64*volume//27+1==1935
assert 1931/(64*volume)==F(30171875,816293376)<F(1,27)
print('PASS: good seed probability at least',lo)
print('PASS: at least 128 light fibres and at most 1931 motif points.')
