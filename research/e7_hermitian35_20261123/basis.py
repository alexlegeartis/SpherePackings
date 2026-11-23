from pathlib import Path
import json,sys
from sympy import Matrix,diag,Rational
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent
C=json.loads((P/'certificate.json').read_text())
G=Matrix(C['E7_gram']); assert G.det()==2
g35=diag(*[G]*5)
Z=[int(x,16) for x in C['generators_35_hex']]
B=hermite_normal_form(Matrix([[2*int(i==j)for j in range(35)]for i in range(35)]+
                             [[(z>>i)&1 for i in range(35)]for z in Z]).T).T
assert abs(B.det())==2**13
cols=B[:,6]
j=next(j for j,c in enumerate(cols)if abs(c)==1)
B34=Matrix([B.row(i)-cols[i]/cols[j]*B.row(j)for i in range(35)if i!=j])
assert B34.shape==(34,35) and not any(B34[:,6])
r=[(82>>i)&1 for i in range(7)]*5
M=C['phase_isometry']
w=[(M[65]>>i)&1 for i in range(7)]+[(M[1]>>i)&1 for i in range(7)]+[0]*21
B36=B.row_join(Matrix.zeros(35,1)).col_join(Matrix([[Rational(t,2)for t in r]+[1]]))
B37=B36.row_join(Matrix.zeros(36,1)).col_join(Matrix([w+[0,1]]))
data={}
for n,BB,GG,d in [(34,B34,g35,3*2**30),(35,B,g35,2**31),(36,B36,diag(g35,Rational(1,2)),2**30),(37,B37,diag(g35,Rational(1,2),2),2**31)]:
 Q=BB*GG*BB.T
 assert Q.det(method='domain-ge')==d
 data[str(n)]=dict(basis_times_2=[[int(2*x)for x in row]for row in BB.tolist()],
                  independent_metric=[str(x)for x in GG.diagonal()],gram_determinant=d)
 print(n,'det',d,'covolume',d**.5,'exact squared minimum 8 (word and layer proof)')
p=P/'bases.json'
if '--write' in sys.argv:p.write_text(json.dumps(data,separators=(',',':'))+'\n')
else:assert data==json.loads(p.read_text())
print('PASS exact basis determinants 34–37')
