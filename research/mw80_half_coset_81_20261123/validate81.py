"""Independent exact algebraic validation for the 81D half-coset proof."""
import json
import sys
from pathlib import Path
from fractions import Fraction
from math import isqrt
import numpy as np
from flint import fmpz_mat as Z, fmpq_mat as Q

root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
if (root/'geometry.npz').exists():
    x = np.load(root/'geometry.npz')
    y = np.load(root/'split_V2^21.npz')
    G=Z(x['G'].tolist()); T=Z(y['T'].tolist())
    A=(Z(y['B'].tolist()), Z(y['C'].tolist()))
else:
    h=json.loads((root/'split_matrices.json').read_text())
    G,T,A=Z(h['G']),Z(h['T']),(Z(h['B']),Z(h['C']))
eye = Z([[int(i==j) for j in range(80)] for i in range(80)])
assert T*T==eye and T*G*T.transpose()==G
assert G.det()==41**12 and not any((A[0]*G*A[1].transpose()).entries())

def rank2(rows):
    ref={}
    for r in rows:
        x=sum((y%2)<<i for i,y in enumerate(r))
        while x:
            p=x.bit_length()-1
            if p in ref: x ^= ref[p]
            else: ref[p]=x; break
    return len(ref)

p=(21382<<64)|1647989595034592450
r=Z([[(p>>j)&1 for j in range(80)]])
assert rank2(A[0].tolist()+A[1].tolist())==42
assert rank2(A[0].tolist()+r.tolist())==rank2(A[1].tolist()+r.tolist())==40
assert not any(x%2 for x in (r*T-r).entries()) and (r*G*r.transpose())[0,0] == 1012
w=Z([list(map(int,(root/'small_rep.txt').read_text().split()))])
assert all(x%2==0 for x in (w-r).entries()) and (w*G*w.transpose())[0,0]==40
rr=(r*G).tolist()[0]
C=Z([[4*G[i,j] for j in range(80)]+[2*rr[i]] for i in range(80)]+[[2*x for x in rr]+[1028]])
assert C.det()==4**82*41**12

def read(name):
    f=list(map(int,(root/name).read_text().split()))
    assert len(f)==40**2+40*80+40+40**2+40
    pos=0
    def take(r,c):
        nonlocal pos
        z=[f[pos+i*c:pos+(i+1)*c] for i in range(r)]; pos+=r*c
        return z
    H=Z(take(40,40));B=Z(take(40,80));D=take(1,40)[0];S=take(40,40);W=take(1,40)[0]
    return H,B,D,S,W

for i in range(2):
    a=A[i]; ht=a*G*a.transpose();h=Q(ht);B=Q(a)
    F=Q(eye+T if i==0 else eye-T)
    assert rank2(a.tolist())==40 and a*T==(a if i==0 else -a)
    L=F*Q(G)*B.transpose()*h.inv()
    assert all(v.denominator==1 for v in L.entries()) and L*B==F
    for name in [f'fixed{i}.in',f'jointK{i}.in']:
        ht,b,d,mu,whole=read(name)
        assert b*G*b.transpose()==ht
        inv=Q(ht).inv()
        assert any(x for x in b.entries()) and rank2(b.tolist()+a.tolist())==40
        if name.startswith('fixed'):assert b==a
        else:
            t=Q(b)*Q(G)*B.transpose()*h.inv()
            assert t*B==Q(b) and all(y.denominator==1 for y in t.entries()) and abs(t.det())==4
            assert rank2(b.tolist())==38 and rank2(b.tolist()+A[1-i].tolist())==40
            assert all(v%2==0 for v in ht.entries()) and all(ht[a,a]%4==0 for a in range(40))
        D=[];U=[[Fraction(int(j==k)) for k in range(40)] for j in range(40)]
        S=1<<20
        for j in range(40):
            dv=Fraction(int(ht[j,j]))-sum(U[j][l]**2*D[l] for l in range(j))
            assert dv>0;D.append(dv)
            assert (dv*S).__floor__()==d[j]
            assert isqrt(int((40*inv[j,j]).floor()))==whole[j]
            for k in range(j+1,40):
                U[k][j]=(Fraction(int(ht[k,j]))-sum(U[j][l]*U[k][l]*D[l] for l in range(j)))/dv
        assert all((U[j][k]*S).__floor__()==mu[j][k] for j in range(40) for k in range(40))
print('PASS exact determinant, eigensublattice saturation and shared-parity index-four embeddings; directed LDL rounding/box verified')
