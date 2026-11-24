"""Check the six missing index-47^2 witnesses in the established cusp basis."""
import json, runpy, sys
from pathlib import Path
here = Path(__file__).resolve().parent
base = (Path(sys.argv[1]) if len(sys.argv) > 1 else
        here.parent.parent / 'scoped_density_round17_20261124' / '92-MW-cusp-jets' / 'verify.py')
g = runpy.run_path(str(base))
rows = json.loads((here / 'witnesses.json').read_text())
assert len(rows) == 6 and {tuple(r[:3]) for r in rows} == {k for k, _ in g['orbit_reps']}
for a,b,c,d,v in rows:
    v = g['np'].array(v,dtype='i8')
    assert v.shape == (92,) and d in (16,18) and g['q'](v) == d
    assert g['case'](g['mm'](v,g['J']) % 47) == (a,b,c)
    print((a,b,c,'norm',d))
assert 18**46 < 16**46 * 47**2
p=g['p']
np=g['np']
mm=g['mm']
T,I,H,A=g['T'],g['I'],g['H'],g['A']
group=[g[k] for k in ('U','C','F')]
for X,k in zip(group,[5,1,1]):
    P=I.copy()
    for _ in range(k):
        P=mm(P,T)
    assert np.array_equal(mm(X,T),mm(P,X))
def leading(v):
    assert g['q'](v)==20
    for j in range(2):
        v=mm(v,H)
        assert not np.any(v%p)
        v=v//p
    a,b=map(int,mm(v,A)%p)
    assert a or b
    return b*pow(a,-1,p)%p if a else p
# Propagate only newly found leading directions. As each isometry
# normalizes T, it permutes K2/K3 and these projective directions.
rows=[g['D']['depth2'],json.loads((here/'norm20.json').read_text())]
directions=[]
for row in rows:
    start=np.array(row,dtype='i8')
    seen={leading(start)}
    work=[start]
    while work:
        v=work.pop()
        for X in group:
            w=mm(v,X)
            i=leading(w)
            if i not in seen:
                seen.add(i)
                work.append(w)
    directions.append(seen)
assert directions==[{0,4,12,35,43,47},{1,10,14,33,37,46}]
assert 20**46<16**46*47**3
assert 48 + len(set.union(*directions))*47**2==26556
print('Index-47^3 lower bound on number of excluded profiles:',26556)
