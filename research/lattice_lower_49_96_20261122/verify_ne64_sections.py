from pathlib import Path
import re, sys
from sympy import Matrix

HERE=Path(__file__).resolve().parent
default=HERE.parent/'novel_routes_20261113/higher_parent/CQ64b.html'
p=Path(sys.argv[1]) if len(sys.argv)>1 else default
t=p.read_text().split('NAME="GRAM"',1)[1].split('NAME="DIVISORS"',1)[0]
data=[int(v) for v in re.findall(r'-?\d+',re.sub('<[^>]*>',' ',t))]
if len(data)==2082:
    assert data[:2]==[64,0]
    G=Matrix.zeros(64)
    for i in range(64):
        for j in range(i+1): G[i,j]=G[j,i]=data[2+i*(i+1)//2+j]
else:
    assert data[0]==64 and len(data)==4097
    G=Matrix(64,64,data[1:])
assert G==G.T and G.det(method="domain-ge")==3**32
for n,K,d in [
 (59,[39,47,48,49,55],296483230216294560),
 (60,[39,47,55,56],118593292086517824),
 (61,[39,47,55],59296646043258912),
 (62,[39,47],22236242266222092),
 (63,[39],4*3**32),
]:
 keep=[i for i in range(64) if i not in K]
 V=G.extract(keep,keep)
 assert V.shape==(n,n) and V[0,0]==12 and V.det(method="domain-ge")==d
 print('PASS Ne64 bare subset',n,d)
