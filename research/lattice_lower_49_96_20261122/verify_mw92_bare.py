"""Exact primitive-kernel proof for MW92 lattice sections in 86–91D."""
from pathlib import Path
from hashlib import sha1
from fractions import Fraction as Q
import json, sys
import numpy as np
from sympy import Matrix as Z, eye
from sympy.matrices.normalforms import hermite_normal_form
from verify_catalog import ball_squared, pi_floor, float_floor

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
S=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'research/packing_release8590_20261117/packing_resume_20261117'
def pin(name, digest):
    p=S/name; d=p.read_bytes()
    assert sha1(f'blob {len(d)}\0'.encode()+d).hexdigest()==digest, name
    return p
z=np.load(pin('parent_certificate.npz','56e2a6f803556a84009b517956960080710610e9'),allow_pickle=False)
G=Z(z['G'].tolist())
assert G.shape==(92,92) and G==G.transpose() and G.det(method="domain-ge")==47**14
invg=G.inv(method="DM")
witnesses=json.loads((HERE/'section_witnesses.json').read_text())
assert sorted(x['dimension'] for x in witnesses)==list(range(86,92))
rows=[]
pi=pi_floor()
for entry in sorted(witnesses,key=lambda x:x['dimension']):
    n=entry['dimension']; k=92-n
    expected=int(entry['det'])
    B=Z(entry['normal_labels'])
    assert B.shape==(k,92)
    # HNF(A)=I is exact surjectivity Z^92 → Z^k. The coarea determinant
    # identity then gives the primitive kernel Gram determinant.
    assert hermite_normal_form(B)==eye(k)
    K=B*invg*B.T
    assert (47*K).tolist()==entry["numerator_gram"]
    assert K.det(method="domain-ge")*47**14==expected
    # Parent minimum 16 is inherited from the published Mordell–Weil theorem.
    atoms=z['section_coordinates']
    arr=np.asarray(B.tolist(),dtype='i8')
    contact=next((list(map(int,w)) for w in atoms if not any(w@arr.T)
                  and int(w@z['G']@w)==16),None)
    if contact is not None:
        v=Z([contact]); assert (v*G*v.transpose())[0,0]==16 and (v*B.transpose()).is_zero_matrix
    D=Q(2**(2*n),expected)
    rows.append({'dimension':n,'kind':'audited lattice section','minimum_squared':16,
        'minimum_attained_by_saved_parent_vector':contact is not None,'minimum_witness':contact,
        'normal_gram':[[str(x) for x in row] for row in K.tolist()],
        'gram_determinant':str(expected),'center_density_squared':str(D),
        'occupied_density_floor':float_floor(D*ball_squared(n,pi),15)})
    print('PASS MW92 section',n,expected,'minimum attained',contact is not None,rows[-1]['occupied_density_floor'])
out=HERE/'mw92_bare_results.json'
if '--emit' in sys.argv:
    out.write_text(json.dumps({'rows':rows},indent=2)+'\n')
else:
    assert json.loads(out.read_text())=={'rows':rows}
# Strict, exact comparison with the preceding unrestricted two-coset result.
old=4*(Q(4337,274)/4)**91/(274*47**13)
new=Q(2**182,182*47**13)
assert new*10_000_000**2>old*10_018_453**2
print('PASS new 91D lattice improves the preceding unrestricted certificate by >0.18453%.')
