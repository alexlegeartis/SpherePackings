from pathlib import Path
import json,sys
from sympy import Matrix,diag
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent
C=json.loads((P/'certificate.json').read_text())
D=json.loads((P/'certificate38.json').read_text())
Z=[int(s,16)for s in C['generators_35_hex']]
B=hermite_normal_form(Matrix([[2*int(i==j)for j in range(35)]for i in range(35)] + [[z>>i&1 for i in range(35)]for z in Z]).T).T
assert abs(B.det())==2**13
bottom=[sum(x,[])+y for x,y in zip(D['horizontal_generators_times_2'],D['vertical_generators_times_2'])]
A=(2*B).row_join(Matrix.zeros(35,3)).col_join(Matrix(bottom))
G=diag(*[Matrix(C['E7_gram'])]*5,1,1,1)
Q=A*G*A.T/4
assert Q.det(method='domain-ge')==2**31
R=dict(basis_coefficients_times_2=A.tolist(),gram_times_4=(4*Q).tolist())
R={k:[[int(x)for x in y]for y in val]for k,val in R.items()}
f=P/'basis38.json'
if '--write' in sys.argv:f.write_text(json.dumps(R,separators=(',',':'))+'\n')
else:assert R==json.loads(f.read_text())
print('PASS exact 38-row lattice basis, determinant 2^31')
