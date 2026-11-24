from joint import jointrows
from second import predicted
from span import admissible
from fractions import Fraction as F
from pathlib import Path
import json
root=Path(__file__).resolve().parent
R,V,M,A,B=jointrows()
assert len(A)==222 and len(B)==441 and len(R)==652
W=lambda f:{int(k):F(v) for k,v in json.loads((root/f).read_text()).items()}
W0=W('anchor-counts.json')
R0=[(j,sum(r[i]*v for i,v in W0.items())-w) for j,(r,w) in enumerate(zip(R,V)) if j<120]
assert all(not z for j,z in R0)
P=predicted(W0)
negative={k:v for k,v in P.items() if v<0}
assert negative=={(4,5):F(-13781219850,7471),(5,4):F(-13781219850,7471),(5,5):F(-50323730375,7471)}
print('Round18 witness has three negative norm12 parity-averaged coefficients:',negative)
T=W('joint-counts.json')
assert len(T)==151 and T[0]==1 and min(T.values())==1
for r,w,m in zip(R,V,M):
 assert sum(r[i]*v for i,v in T.items())==w,m
print('Replacement exact histogram:',len(T),'positive orbit weights;',len(R),'verified rows.')
for j in T:
 if j<len(A):
  assert all(admissible(10,*z) for z in A[j])
 else:
  assert all(admissible(12,*z) for z in B[j-len(A)])
print('Every used orbit additionally passes exhaustive arithmetic sign/shell check through norm16 in its complete span.')
from affine import extra_rows,advance,shift
assert shift(advance((F(10),F(-4),F(2),F(2)),-1),-1,0)==(F(12),F(-4),F(-4),F(-4))
assert (-4,2,2) in A[54] and (-4,-4,-4) in B[306-len(A)]
assert T.get(54,0)==0 and T.get(306,0)==227401200
r2,v2,m2,a2,b2=extra_rows()
assert len(r2)==873 and a2==A and b2==B
T2=W('affine-counts.json')
assert len(T2)==153 and T2[0]==1 and all(x>0 for x in T2.values())
for r,w,m in zip(r2,v2,m2):
 assert sum(r[i]*v for i,v in T2.items())==w,m
for j in T2:
 if j<len(A):
  assert all(admissible(10,*z) for z in A[j])
 else:
  assert all(admissible(12,*z) for z in B[j-len(A)])
print('PASS all 873 rows, including all affine-unit identifications in the two shells.')
from joint3 import jointrows3
r3,v3,m3,a3,b3,c3=jointrows3()
assert a3==A and b3==B and len(c3)==728 and len(r3)==1953
T3=W('joint3-counts.json')
assert {j:z for j,z in T3.items() if j<len(A)+len(B)}==T2
assert len(T3)==298 and all(z>0 for z in T3.values())
for r,w,m in zip(r3,v3,m3):
 assert sum(r[i]*v for i,v in T3.items())==w,m
for j in T3:
 if j>=len(A)+len(B):
  assert all(admissible(14,*z) for z in c3[j-len(A)-len(B)])
print('PASS 1953 marked three-shell identities, and all 221 additional two-shell affine-unit rows.')
