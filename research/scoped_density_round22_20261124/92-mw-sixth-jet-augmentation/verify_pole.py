"""Exact finite-field rational section identity for the new sixth-jet witness."""
import json,runpy
from pathlib import Path
here=Path(__file__).resolve().parent
release=runpy.run_path(str(here/'verify.py'))
base=release['d']['g']; p=base['p']; np=base['np']
d=json.loads((here/'pole.json').read_text())
x,y,a=(d[k]for k in ('xnum','ynum','pole'))
r,minus,powp,mul,gcddeg=(base[k]for k in ('trim','sub','power','mul','gcddeg'))
g=lambda x: int(x)%p+(int(x)//p)*p
assert all(g(c)==c and 0<=c<p*p for c in x+y+[a])
F=[base['amin'](0,a),1]; h=base['f']
assert len(r(x[:]))==27 and len(r(y[:]))==40
assert minus(powp(y,2),powp(x,3))==mul(powp(F,6),h)
def evalp(P):
    out=0
    for c in reversed(P):
        out=base['aplus'](base['times'](out,a),c)
    return out
assert evalp(x) and evalp(y)
assert a==40+34*p and a>=p
w=np.array(release['R'][0]['w'],dtype='i8')
Gram,W,mm=(base[k]for k in ('G','W','mm'))
pairings=mm(mm(w,Gram),W.T)
# All fibres are irreducible. The new section intersects O four
# times at infinity and once at a. The reference norm-16 sections are
# integral and specialize away from O at both sites.
F2,F3=powp(F,2),powp(F,3)
got=[13-gcddeg(minus(x,mul(s,F2)),minus(y,mul(t,F3))) for s,t in base['sections']]
assert np.array_equal(pairings,got)
print('Height pairings identify the exact lattice row with its rational section; unique finite pole is 40+34i.')
