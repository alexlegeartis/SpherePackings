from itertools import product
from math import isqrt
from sympy import Matrix
from fractions import Fraction as F
import numpy as np

G=2*np.eye(6,dtype=np.int64)
for i,j in [(0,1),(1,2),(2,3),(3,4),(2,5)]:G[i,j]=G[j,i]=-1
inv=Matrix(G.tolist()).inv()
v=(2,1,0,2,1,0)
out=[]
for a in (v,tuple((2*t)%3 for t in v)):
    # Q(x,y) >= 3/4 Q(x), 3/4 Q(y), so enumerate Q(x)<=24.
    B=[isqrt(int(24*inv[i,i])) for i in range(6)]
    S=np.array(list(product(*[[j for j in range(-b,b+1) if j%3==a[i]] for i,b in enumerate(B)])),dtype=np.int64)
    n=np.sum((S@G)*S,axis=1)
    T=S[n<=24]; m=n[n<=24]
    assert len(T) and min(m)==12
    pair=m[:,None]+m[None,:]+T@G@T.T
    assert int(pair.min())==18
    print('legal residue',a,'lifts',len(T),'minimum numerator',int(pair.min()))
assert Matrix(G.tolist()).det()==3
# E6 has even squared norms and is positive definite; the two
# pure vectors of a root have squared length 2. For nonzero (x,y),
# when both are nonzero, Q(x,y)>=Q(x)/2+Q(y)/2>=2.
covolume2=9*F(3,4)**6
assert covolume2==F(81,64)**2
assert F(3,64)/F(81,64)==F(1,27)
print('matched center density 1/27')
