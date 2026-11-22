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
phase=json.loads((here/'phase.json').read_text())
raw=phase['raw']
assert len(raw)==72 and set(raw)<={0,1}
assert P['syndrome'](raw)==tuple(phase['syndrome'])
P['prove'](raw,6)
H=np.vstack([P['rawD4'],raw])
for z in product(range(2),repeat=4):
    r=(np.array(z)@H[:4]+H[4])%2
    a,b,t=P['syndrome'](r)
    S=(a,b,a^b)
    assert all(S)
    assert all(t!=x^y for p,q in combinations(S,2)
               for x in P['M4'][p] for y in P['M4'][q])
    # At least two nonzero Leech blocks (distinct syndrome A labels), and an
    # actual raw squared norm of 8 would therefore require the excluded 4+4.
    # Raw norms even imply horizontal squared depth at least 5.

R=Matrix([[1,-1,0,0],[0,1,-1,0],[0,0,1,-1],[0,0,1,1]])
assert abs(R.det())==2 and all(sum(r)%2==0 for r in R.tolist())
assert R*R.T==P['V0']
half=Matrix([Q(1,2),0,0,0]); U=R*half
V= (P['V0'].row_join(U)).col_join(U.T.row_join(Matrix([[2]])))
assert V.is_positive_definite and V.det()==7
assert 2-(half.T*half)[0]==Q(7,4)
# k odd +/-1: x in D4, vertical squared length 2+||x||² +/- x1.
# If x!=0, at least one other coordinate is nonzero (contributes >=1),
# or |x1| is an even positive integer (contributes >=2). The x1 term n(n+/-1)
# is nonnegative for all integers. Thus vertical >=3, which combines
# with horizontal >=5. At x=0 use depth6 and vertical length²2.
# k=+/-2: vertical²=7+||x+/-e1||²>=8 by odd coordinate sum.
# |k|>=3: vertical²>=9*7/4>8. k=0 is the verified 76D D4 extension.

W=2*Matrix(H.tolist())*P['Ti']
assert all(w.q==1 for w in W)
G=P['GG']; X=W*G/2; N=W*G*W.T/4+V
GG=G.row_join(X.T).col_join(X.row_join(N))
assert GG.det(method='domain-ge')==7
num=4*GG
assert all(a.q==1 for a in num)
data={'dimension':77,'minimum_squared':8,'gram_determinant':7,
      'center_density_squared':str(Q(2**77,7)),
      'numerator_denominator':4,
      'gram_numerator':[[int(a) for a in row] for row in num.tolist()]}
if '--emit' in sys.argv:
    (here/'gram77.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
else:
    assert data==json.loads((here/'gram77.json').read_text())
print('PASS: 77D lattice, squared minimum 8, Gram determinant 7, squared center density 2^77/7')
