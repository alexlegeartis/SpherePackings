"""Independent determinant, positive definiteness, and congruence audit."""
from pathlib import Path
from sympy import Matrix
from sympy.polys.matrices import DomainMatrix
import gzip,json,hashlib
P=Path(__file__).resolve().parent
D=json.loads(gzip.decompress((P/'gram96.json.gz').read_bytes()))
H=Matrix(D['horizontal_gram'])
T=Matrix(D['glue_basis_from_horizontal_and_D4'])
G=Matrix(D['gram'])
B=Matrix([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]])
C=Matrix(D['cusp_map']+[[-1,0],[-26,-33],[-36,-35],[-35,-37]])
def det(A): return int(DomainMatrix.from_Matrix(A).det())
assert H.shape==(92,92) and G.shape==T.shape==(96,96)
assert G==T*Matrix.diag(H,B)*T.T
assert H==H.T and G==G.T and all(int(G[i,i])%2==0 for i in range(96))
assert all(int(u)%47==0 for u in T*C)
assert abs(det(T))==47**2 and det(H)==47**14 and det(G)==4*47**18
# Bareiss pivots are the leading principal minors, hence positive implies SPD.
A=G.tolist(); last=1
for k in range(96):
 p=A[k][k];assert p>0
 for i in range(k+1,96):
  for j in range(k+1,96):
   q,r=divmod(A[i][j]*p-A[i][k]*A[k][j],last)
   assert r==0;A[i][j]=q
 last=p
assert last==4*47**18
print('PASS independent congruence, SPD and determinants; Gram SHA256',hashlib.sha256(json.dumps(D['gram'],separators=(',',':')).encode()).hexdigest())
