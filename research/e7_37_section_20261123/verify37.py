from pathlib import Path
import json,sys
from sympy import Matrix,Rational,diag
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent
from hashlib import sha256
from subprocess import run
R=P.parent/'e7_lattice38_20261123'
for name in ('verify35.py','verify38.py','basis.py'):
    run([sys.executable,'-B',str(R/name)],check=True)
C=json.loads((R/'certificate.json').read_text())
A=json.loads((R/'certificate38.json').read_text())
G=Matrix(C['E7_gram']); D=diag(*([G]*5),1,1,1)
Z=[int(x,16) for x in C['generators_35_hex']]
H=hermite_normal_form(Matrix([[2*int(i==j) for j in range(35)] for i in range(35)]+
                            [[(z>>i)&1 for i in range(35)] for z in Z]).T).T
assert H.det()==2**13 and (H*D[:35,:35]*H.T).det(method='domain-ge')==2**31
V=Matrix(A['vertical_generators_times_2'])
F=H.row_join(Matrix.zeros(35,3)).col_join(Matrix([sum(w,[])+v for w,v in zip(A['horizontal_generators_times_2'],V.tolist())])/2)
w=Matrix([0]*35+[Rational(1,2),-Rational(1,2),Rational(1,2)])
pair=F*D*w
assert pair==Matrix([0]*35+[0,1,0])
assert (w.T*D*w)[0]==Rational(3,4) and F.det()==2**13*V.det()/8
s=[i for i in range(38) if i!=36]
F37=F[s,:]
Q37=F37*D*F37.T
assert Q37.det(method='domain-ge')==3*2**29
v2=[4]+[0]*37
assert (Matrix([v2])*D*Matrix(v2))[0]==32
# Doubled first simple root belongs to horizontal period, with exact integral row coefficients.
assert all(z.q==1 for z in Matrix([[2]+[0]*34])*H.inv())
assert Rational(2)**37/Q37.det(method="domain-ge")==Rational(256,3)
out=dict(
    dimension=37,
    basis_times_2=[[int(2*x) for x in r] for r in F37.tolist()],
    gram_times_4=[[int(4*x) for x in r] for r in Q37.tolist()],
    norm2=8,
    gram_determinant=3*2**29,
    center_density_squared='256/3',
    parent_gram_determinant=2**31,
    normal_vertical_times_2=[1,-1,1],
    parent_normal_pairings=[0]*35+[0,1,0],
)
p=P/'section37.json'
if '--write' in sys.argv:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
else:assert json.loads(p.read_text())==out
print('PASS 37D exact basis, determinant 3*2^29 and squared minimum 8 inherited from certified parent; center density 16/sqrt3')
