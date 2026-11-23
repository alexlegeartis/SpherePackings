import json,math,itertools
from sympy import Matrix
from pathlib import Path
G=Matrix(json.loads(Path('certificate.json').read_text())['E7_gram'])
L,D=G.LDLdecomposition(hermitian=False)
S=math.lcm(*[int(x.q) for x in L])
T=math.lcm(*[int(x.q) for x in D.diagonal()])
assert all((x>0)for x in D.diagonal())
print(S,T)
print(*[int(x*T)for x in D.diagonal()])
for i in range(7):print(*[int(y*S)for y in L[i,:]])
