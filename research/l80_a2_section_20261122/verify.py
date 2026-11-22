#!/usr/bin/env python3
"""Exact 82D restriction of the certified 83D L80 + A3 lattice."""
from fractions import Fraction as Q
from functools import reduce
from hashlib import sha1
import json
from math import factorial,lcm
from pathlib import Path
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'packing83_85_deep_verified_20261117/deep_extension/verified83.json'
SOURCE_BLOB='a8a4ff05846d126e29b0580967989832dfea4b00'
def bareiss(A):
    B=[row[:] for row in A]
    d=1
    for j in range(len(B)-1):
        k=B[j][j];assert k>0
        for i in range(j+1,len(B)):
            for c in range(j+1,len(B)):
                num=k*B[i][c]-B[i][j]*B[j][c]
                q,r=divmod(num,d);assert not r
                B[i][c]=q
            B[i][j]=0
        d=k
    assert B[-1][-1]>0
    return B[-1][-1]
def det(A):
    assert all(len(r)==len(A) for r in A)
    d=lcm(*(x.denominator for r in A for x in r))
    return Q(bareiss([[int(x*d) for x in r]for r in A]),d**len(A))
def arctan_pair(x, n):
    # n terms start with a positive term. Even n => last term negative,
    # valid lower bound, and the following positive term an upper bound.
    assert n%2==0 and 0<x<1
    y=sum((-1)**i*x**(2*i+1)/(2*i+1) for i in range(n))
    return y,y+x**(2*n+1)/(2*n+1)
def verify():
    assert __debug__
    raw=SOURCE.read_bytes()
    assert sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==SOURCE_BLOB
    src=json.loads(raw)
    assert src['dimension']==83 and src['determinant']==4 and src['minimum_squared_at_least']==8
    A=[[Q(x)for x in row]for row in src['gram']]
    assert len(A)==83 and all(A[i][j]==A[j][i]for i in range(83)for j in range(83))
    assert det(A)==4 and A[3][3]==8
    # Removing the final basis vector gives an actual sublattice, not a
    # projection. Its shortest length is bounded below by the 83D parent
    # minimum; the retained fourth vector attains it.
    prefix=lambda n:[r[:n]for r in A[:n]]
    assert det(prefix(80))==1 and det(prefix(81))==2 and det(prefix(82))==3
    m=82
    center=Q(2**m,3)
    old=Q(2**86,55)
    assert center/old==Q(55,48)
    a,_=arctan_pair(Q(1,5),70)
    _,b=arctan_pair(Q(1,239),24)
    pi=16*a-4*b
    lower_square=pi**m*center/factorial(m//2)**2
    floor='9.17033181889404e-18'
    assert Q(floor)**2<lower_square
    gain='7.0436048222094'
    assert (1+Q(gain)/100)**2<center/old
    return dict(dimension=m,kind='explicit lattice',minimum_squared=8,
                gram_determinant=3,previous_center_density_squared=str(old),
                center_density_squared=str(center),ratio_squared='55/48',
                occupied_density_floor=floor,gain_over_repository_percent_floor=gain,
                parent_83_gram_blob=SOURCE_BLOB)
if __name__=='__main__':
    result=verify()
    assert result==json.loads((HERE/'results.json').read_text())
    print(json.dumps(result,indent=2))
