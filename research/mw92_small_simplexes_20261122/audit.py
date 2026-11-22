"""Independent exact matrix inverse and separation/covolume audit."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from fractions import Fraction as F
import json,numpy as np
from sympy import Matrix,eye
from sympy.polys.matrices import DomainMatrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent
q=next(a for a in (P.parent/'mw92_sections_20261115'/'parent_certificate.npz',
                  P.parent/'lattice_catalog85_96'/'parent0.npz') if a.is_file())
assert sha256(q.read_bytes()).hexdigest()=='3c2a8629853cbd3f4ac37ad128a04b8a4c21db67a41888fed766ae840e0013f0'
G=Matrix(np.load(q,allow_pickle=False)['G'].tolist())
Gi=DomainMatrix.from_Matrix(G).to_field().inv().to_Matrix()
assert G*Gi==eye(92) and G.det(method='domain-ge')==47**14
rows=json.loads((P/'verified_results.json').read_text())
for j in rows:
    n=j['dimension'];w=json.loads((P/f'band_winner{n}.json').read_text())
    k=92-n
    B=Matrix(w['B']);assert B.shape==(k,92) and hermite_normal_form(B)==eye(k)
    E=Gi*B.T;H=B*E
    assert all(H[:i,:i].det()>0 for i in range(1,k+1))
    assert H*47==Matrix(w['K'])
    det=47**14*H.det()
    assert det==F(j['period_determinant'])
    labels=[tuple(x)for x in w['L']]
    assert len(labels)==len(set(labels))==j['translates']
    C=E*H.inv()
    loss=[]
    for a,b in combinations(labels,2):
        z=Matrix(a)-Matrix(b)
        v=C*z
        assert B*v==z
        loss.append((v.T*G*v)[0,0])
    beta=max(loss)
    assert beta==F(j['beta']) and 0<beta<16
    sq=F(((16-beta)/4)**n)*len(labels)**2/F(det)
    assert sq==F(j['center_density_squared'])
    print(n,'PASS exact normal representatives, covolume and density',det,beta)
