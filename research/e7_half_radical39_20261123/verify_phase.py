"""Verify the local enumeration bounds and replay all horizontal phase depths."""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from tempfile import TemporaryDirectory
from subprocess import run
import json
P=Path(__file__).resolve().parent
A=P.parent/'e7_hermitian35_20261123'/'certificate.json'
d=A.read_bytes()
assert sha256(d).hexdigest()=='480e7bba7424119eca9eef0339b7b4b3a6c93b7c3ee7f8db575f867634d0c497'
q=json.loads(d)
G=q['E7_gram']
H=[list(map(F,row))+[F(i==j)for j in range(7)]for i,row in enumerate(G)]
for i in range(7):
 a=H[i][i];assert a>0;H[i]=[x/a for x in H[i]]
 for j in range(7):
  if j==i:continue
  a=H[j][i]
  H[j]=[x-a*y for x,y in zip(H[j],H[i])]
assert [2*H[i][i+7]for i in range(7)]==[24,7,12,4,15,8,3]
b=[12,7,9,5,10,7,4]
assert all(b[i]**2<=14*H[i][i+7]<(b[i]+1)**2 for i in range(7))
assert q['field_polynomial']==67
with TemporaryDirectory() as t:
 exe=Path(t)/'phase_depth'
 run(['g++','-O3','-std=c++17',str(P/'phase_depth.cpp'),'-o',str(exe)],check=True)
 run([str(exe)],input=' '.join(map(str,q['phase_isometry']+q['alpha']))+'\n',text=True,check=True)
