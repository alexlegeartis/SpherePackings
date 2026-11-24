"""Integer checks for the 12D anisotropic E6 x E6 quadratic-label obstruction."""
from collections import Counter
from itertools import product
from math import isqrt
from fractions import Fraction as F
import numpy as np

G=2*np.eye(6,dtype=np.int64)
for i,j in [(0,1),(1,2),(2,3),(3,4),(2,5)]:
    G[i,j]=G[j,i]=-1
I=np.array([[4,5,6,4,2,3],[5,10,12,8,4,6],[6,12,18,12,6,9],
    [4,8,12,10,5,6],[2,4,6,5,4,3],[3,6,9,6,3,6]],dtype=np.int64)
assert np.array_equal(G@I,3*np.eye(6,dtype=np.int64))
pw=3**np.arange(6); V=(np.arange(729)[:,None]//pw)%3
rad=np.array([2,1,0,2,1,0])
assert set(np.flatnonzero(((V@G)%3==0).all(axis=1)))=={0,rad@pw,((2*rad)%3)@pw}
# Find literal lift witnesses up to squared norm 8, using inverse Gram.
B=[isqrt(8*int(v)//3) for v in np.diag(I)]
X=np.array(list(product(*[range(-k,k+1) for k in B])),dtype=np.int64)
norms=np.sum((X@G)*X,axis=1)
ids=(X%3)@pw
d=np.full(729,99,dtype=np.int64)
np.minimum.at(d,ids,norms)
assert set(np.flatnonzero(d>=10))=={rad@pw,((2*rad)%3)@pw}
quotient_types=(np.sum((V@G)*V,axis=1))%3
shifts=((V[:,None,:]+np.arange(3)[None,:,None]*rad)%3)@pw
assert np.all(d[shifts[quotient_types!=0]].min(axis=1)<=4)
roots=X[norms==2]
assert len(roots)==72
def rank3(A):
    a=np.array(A,dtype=np.int64)%3
    r=0
    for j in range(a.shape[1]):
        cand=np.flatnonzero(a[r:,j])
        if len(cand)==0:continue
        i=r+int(cand[0]);a[[r,i]]=a[[i,r]]
        a[r]=(a[r,j]*a[r])%3
        a[r+1:]=(a[r+1:]-a[r+1:,j,None]*a[r])%3
        r+=1
    return r
def ev(R):
    return [[r[i]*r[j] for i in range(len(r)) for j in range(i,len(r))] for r in R]
assert rank3(ev(roots))==21
A4roots=roots[((roots[:,4:]%3)==0).all(axis=1),:4]
assert len(A4roots)==20 and rank3(ev(A4roots))==10
# Zero-sum five-coordinate model of the first four simple roots.
# A literal short lift is obtained by centering and, if needed,
# changing one ±1 by ∓3 to remove the coordinate sum ±3.
allowed={}
V4=[v for v in product(range(3),repeat=5) if sum(v)%3==0]
assert len(V4)==81
for v in V4:
    cen=[x if x!=2 else -1 for x in v]
    s=sum(cen); depth=sum(x*x for x in cen)+abs(s)
    if s:
        assert abs(s)==3 and s//3 in cen
        cen[cen.index(s//3)]-=s
    assert sum(cen)==0 and sum(x*x for x in cen)==depth
    if depth<6:continue
    k=(v.count(0),min(v.count(1),v.count(2)),max(v.count(1),v.count(2)))
    allowed.setdefault(k,[]).append(v)
A,B=(allowed[k] for k in [(2,0,3),(0,1,4)])
assert set(allowed)=={(2,0,3),(0,1,4)} and (len(A),len(B))==(20,10)
def character(u,S):
    return sum(2-3*(sum(i*j for i,j in zip(u,s))%3!=0) for s in S)
z=[character((1,0,0,0,0),S) for S in [A,B]]
w=[character((1,1,0,0,0),S) for S in [A,B]]
assert z==[4,-10] and w==[-14,2]
assert [-(15*x+20*y)//11 for x,y in zip(z,w)]==[20,10]
assert 1+2*F(35,11)==F(81,11) and F(81,11)<8
# Mass test on an 81-point coset: degree-two nonzero => >=27 bad
# pairs, each holding at most 3 sites. Remaining pairs hold at most 6.
assert (27*3+54*6)//2==202
# All projected labels lying in one affine isotropic plane:
# 9 labels, 3 possible horizontal residues and 7 per vertical code.
assert 9*3*7==189<203
assert 9*202==1818 and 1935-1818==117
assert 3*81//2==121 and 9*121==1089
print('E6 lift gates, A4 root span, A4 exact Fourier dual, and occupancy thresholds passed.')
