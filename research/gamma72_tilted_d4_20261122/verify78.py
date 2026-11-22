"""Replay: python verify.py [path/to/the/frozen/gamma72_depth6_20261122]."""
from pathlib import Path
from itertools import combinations, product
from hashlib import sha256
from sympy import Matrix,Rational as Q
import runpy, sys, json, numpy as np

here=Path(__file__).resolve().parent
base=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here.parent/'gamma72_depth6_20261122'
pins={
    'verify.py':'ca65ade908d1399acb93c549f304dcf565f731d347a9ebe525cf1e5a28c7d2c6',
    'certificate.json':'60d65708796ece1a76afd974d77e509fbe01dd2ceef37fbffbf6dd7dfffd8870',
    'numerators.json':'c6a30964135f477fc361bd2c9e6350a4041524d763db2fe901fc2e8abde5fd59',
    'results.json':'46393d59baa87284c98bb4f45920f401f78a0ad66c67c6e612739d397654ba00',
}
for f,h in pins.items():
    assert sha256((base/f).read_bytes()).hexdigest()==h, f
P=runpy.run_path(str(base/'verify.py'))
phases=json.loads((here/'phases.json').read_text())
assert len(phases)==2
for phase in phases:
    raw=phase['raw']
    assert len(raw)==72 and set(raw)<={0,1}
    assert P['syndrome'](raw)==tuple(phase['syndrome'])
    P['prove'](raw,6)
H=np.vstack([P['rawD4']]+[p['raw'] for p in phases])
for z in product(range(2),repeat=6):
    if z[-2:]==(0,0): continue
    r=(np.array(z)@H)%2
    a,b,t=P['syndrome'](r)
    S=(a,b,a^b)
    assert all(S)
    assert all(t!=x^y for p,q in combinations(S,2)
               for x in P['M4'][p] for y in P['M4'][q])

R=Matrix([[1,-1,0,0],[0,1,-1,0],[0,0,1,-1],[0,0,1,1]])
assert abs(R.det())==2 and all(sum(r)%2==0 for r in R.tolist())
assert R*R.T==P['V0']
tilt=Matrix([[Q(1,2),Q(1,2)],[0,0],[0,0],[0,0]])
U=R*tilt
D=Matrix([[Q(7,4),Q(-3,4)],[Q(-3,4),Q(7,4)]])
V= (P['V0'].row_join(U)).col_join(U.T.row_join(D+tilt.T*tilt))
assert V.is_positive_definite and V.det()==10
# Residual vertical height squared is m²+n²+3/4(m-n)². If |m| or |n|
# exceeds 2, it is >=9. Enumerate remaining pairs below height squared 8:
short={}
for m,n in product(range(-2,3),repeat=2):
    if m==n==0: continue
    h=(Matrix([[m,n]])*D*Matrix([m,n]))[0]
    assert h==m*m+n*n+Q(3,4)*(m-n)**2
    if h>=8: continue
    # Exact distance² of (m+n)e1/2 from D4, using coordinate parity.
    t=(m+n)%4
    d=Q(0) if t==0 else Q(1) if t==2 else Q(1,4)
    q=h+d
    if (m%2,n%2)==(0,0):
        assert q>=8
    elif q<3:
        assert h==Q(7,4) and sorted([abs(m),abs(n)])==[0,1]
        # A singleton layer has squared vertical length 2+||x||²±x1,
        # only x=0 is below 3, and uses the unshifted depth-6 phase.
    else:
        assert q>=3
# For each nonzero parity phase, all D4 shifts have depth at least 5.
# The m=n=0 layer is the verified old 76D lattice.
W=2*Matrix(H.tolist())*P['Ti']
assert all(w.q==1 for w in W)
G=P['GG']; X=W*G/2; N=W*G*W.T/4+V
GG=G.row_join(X.T).col_join(X.row_join(N))
assert GG.det(method='domain-ge')==10
num=4*GG
assert all(a.q==1 for a in num)
data={'dimension':78,'minimum_squared':8,'gram_determinant':10,
      'center_density_squared':str(Q(2**78,10)),
      'numerator_denominator':4,
      'gram_numerator':[[int(a) for a in row] for row in num.tolist()]}
if '--emit' in sys.argv:
    (here/'gram78.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
else:
    assert data==json.loads((here/'gram78.json').read_text())
print('PASS: 78D lattice, squared minimum 8, Gram determinant 10, squared center density 2^78/10')
