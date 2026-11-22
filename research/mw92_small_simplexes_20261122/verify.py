from pathlib import Path
from fractions import Fraction as Q
from math import isqrt, factorial
from decimal import Decimal as D, localcontext
from hashlib import sha256
from itertools import combinations
import numpy as np,json
from sympy import Matrix as M
from sympy.matrices.normalforms import hermite_normal_form as hnf
P=Path(__file__).resolve().parent
par=next(path for path in (P/'parent0.npz',P.parent/'mw92_sections_20261115/parent_certificate.npz',P.parent/'lattice_catalog85_96/parent0.npz') if path.is_file())
assert sha256(par.read_bytes()).hexdigest()=='3c2a8629853cbd3f4ac37ad128a04b8a4c21db67a41888fed766ae840e0013f0'
G=np.load(par,allow_pickle=False)['G']; n0=92
assert G.shape==(n0,n0) and G.dtype.kind=='i'
F=np.rint(47*np.linalg.inv(G)).astype('i8')
assert 92*int(abs(F).max())*int(abs(G).max())<2**63 and np.array_equal(G@F,47*np.eye(92,dtype='i8'))
F=M(F.tolist())
def atan(b,n):
    return sum((Q((-1)**i,(2*i+1)*b**(2*i+1)) for i in range(n)),Q(0))
pi=16*atan(5,60)-4*atan(239,61)
def rstr(x,digits=17):
    with localcontext() as ctx:
        ctx.prec=55
        z=(D(x.numerator)/D(x.denominator)).sqrt()
        q=Q(10)**(z.adjusted()+1-digits)
        t=x/q**2; k=isqrt(t.numerator//t.denominator)
        assert (k*q)**2<=x<((k+1)*q)**2
        return f'{D(k*q.numerator)/D(q.denominator):.{digits-1}E}'
old={int(k):v for k,v in json.loads((P/'baseline.json').read_text()).items()}
def check(path):
    j=json.loads(path.read_text());n=j['n'];k=92-n
    B=M(j['B']);S=B*F*B.T
    assert B.shape==(k,92) and hnf(B)==M.eye(k)
    assert all(S[:i,:i].det()>0 for i in range(1,k+1))
    assert S.tolist()==j['K']
    assert len(j['L'])==len(set(map(tuple,j['L']))) and all(len(l)==k for l in j['L'])
    V=[M([x]) for x in j['L']];Mi=len(V)
    T=47*S.inv();beta=max(Q(str(((a-b)*T*(a-b).T)[0,0]))for a,b in combinations(V,2))
    det=int(S.det())*47**(14-k)
    assert 0<=beta<16
    dc=Mi*Mi*((16-beta)/4)**n/det
    assert type(dc)==Q
    r=old[n];prev=Q(r['M']**2,int(r['period_determinant']))*((16-Q(r['beta']))/4)**n
    assert dc>prev
    gt=dc/prev
    gain=isqrt((10**12*gt.numerator)//gt.denominator)-10**6
    m=n//2
    vol=pi**m/factorial(m) if n%2==0 else 2*4**m*pi**m*factorial(m)/factorial(n)
    d=rstr(dc*vol**2)
    print(n, 'translates',Mi, 'det',det,'beta',beta,'occupied lower',d,'gain ≥',gain/10000,'%')
    return dict(dimension=n, translates=Mi, period_determinant=str(det), beta=str(beta),center_density_squared=str(dc),occupied_density_floor=d,gain_percent_floor=f'{gain//10000}.{gain%10000:04d}')
if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser(description="Exact rational replay for the MW92 periodic sections.")
    a.add_argument("--write",action="store_true")
    args=a.parse_args()
    paths=sorted(P.glob('band_winner*.json'))
    j=[check(f)for f in paths]
    assert [r['dimension'] for r in j] in ([86,87],[86,87,88],[86,87,89],[86,87,88,89])
    f=P/'verified_results.json'
    if args.write:f.write_text(json.dumps(j,indent=2)+'\n')
    else:assert j==json.loads(f.read_text())
