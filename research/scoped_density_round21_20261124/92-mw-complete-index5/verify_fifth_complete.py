"""Replay exact completion of the MW92 index-p^5 obstruction."""
import json,runpy,sys
from pathlib import Path
here=Path(__file__).resolve().parent
before=Path(sys.argv.pop(1)) if len(sys.argv)>1 else (
    here.parent.parent/'scoped_density_round20_20261124'/'92-mw-cusp-index4-5'/'verify_fifth.py'
)
d=runpy.run_path(str(before))
p,np,J5,g,mm,encode,dec,group=(d[k] for k in ('p','np','J3','g','mm','encode','dec','group'))
left=set(map(int,d['leftover']))
assert len(left)==24 and sum(d['counts'])==234171672
rows=json.loads((here/'terminal.json').read_text())
assert len(rows)==2
expected=[(2,27,21),(6,24,0)]
for row,label in zip(rows,expected):
    assert tuple(row['case'])==label
    w=np.array(row['w'],dtype='i8')
    assert g['q'](w)==row['norm']==24
    z=mm(w,J5)%p
    assert np.array_equal(z,row['jet']) and not np.any(z[:4]) and np.any(z[4:6])
    key=int(encode(z[4:].reshape(1,6))[0])
    assert key==sum(c*p**(2-j) for j,c in enumerate(label))
    seen={key}; todo=[key]
    while todo:
        x=dec(np.array([todo.pop()],dtype='i8'),3)
        for M in group:
            y=mm(x,M)%p
            assert not np.any(y[:,:4]) and np.any(y[:,4:6])
            v=int(encode(y[:,4:])[0])
            if v not in seen:
                seen.add(v);todo.append(v)
    assert len(seen)==12 and seen<=left
    left-=seen
    print('Depth-two norm-24 terminal witness',label,'closes',len(seen)*p**2,'primitive lifts.')
assert not left
mixed=d['mixed']; total=d['total']
assert sum(d['counts'])+sum(12*p*p for _ in rows)+mixed==total==234330768
assert 24**46<16**46*47**5
print('ALL',total,'translation-invariant index-47^5 sublattices contain a nonzero vector of squared norm <=24.')
