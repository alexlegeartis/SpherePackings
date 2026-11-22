"""Verify the bare 79D MW80 section using a pinned existing parent and its label."""
from pathlib import Path
from hashlib import sha1
from fractions import Fraction as Q
from math import gcd
import json, sys
import numpy as np
from sympy import Matrix
from verify_catalog import ball_squared, pi_floor

HERE=Path(__file__).resolve().parent
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parent/'packing78_verified_20261117'/'mw80_normal162_20261117'
def pin(name,h):
    p=BASE/name
    if not p.exists(): p=BASE/Path(name).name
    a=p.read_bytes()
    assert sha1(f'blob {len(a)}\0'.encode()+a).hexdigest()==h, name
    return p
P=np.load(pin('parent/parent80.npz','1786a085877669bd56b9f0bcb078c658ceb95358'),allow_pickle=False)
G=Matrix(P['G'].tolist()); assert G.shape==(80,80) and G==G.T and G.det(method='domain-ge')==41**12
inv=G.inv(method='DM')
raw=json.loads(pin('candidate78.json','14e0068d01044d31486f382cd2c4e3aca7354b94').read_text())
out=json.loads((HERE/'mw80_79.json').read_text())
assert [x['normal_row'] for x in out]==raw['normal_rows']
d=162*41**11
for rec in out:
    b=rec['normal_row']; B=Matrix([b])
    assert gcd(*b)==1
    assert (B*inv*B.T)[0,0]*G.det(method='domain-ge')==d
    assert int(rec['section_det'])==d
    v=Matrix([rec['norm14_witness']]); assert (v*G*v.T)[0,0]==14 and (v*B.T)[0,0]==0
    q=Q(7,2)**79/d
    assert q==Q(rec['center_density_squared'])
    assert Q(rec['occupied_floor'])**2 < q*ball_squared(79,pi_floor())
    assert q > Q(2**79,1)/Q(27,4)
print('PASS both primitive 79D MW80 sections: determinant 162·41^11 and attained squared minimum 14.')
