"""Exact augmented MW92 sixth-jet short-witness screen."""
import json,runpy,sys
from pathlib import Path
here=Path(__file__).resolve().parent
before=Path(sys.argv.pop(1)) if len(sys.argv)>1 else (
    here.parent.parent/'scoped_density_round21_20261124'/'92-mw-cusp-index6'/'verify.py'
)
d=runpy.run_path(str(before))
p,np,J6,g,mm,encode,dec,group=(d[k]for k in ('p','np','J3','g','mm','encode','dec','group'))
R=json.loads((here/'additional.json').read_text())
rows=[]
for row in R:
    wi=np.array(row['w'],dtype='i8')
    assert g['q'](wi)==row['norm'] and 16<=row['norm']<=26
    y=mm(wi,J6)%p
    assert np.array_equal(y,row['jet'])
    dep=row['depth']
    assert 2<=dep<6 and not np.any(y[:2*dep]) and np.any(y[2*dep:2*dep+2])
    rows.append(y)
X=np.vstack((d['X'],rows))
active={}
def allowed(z,L):
    a=z//p**L;r=z%p**L
    ans=np.ones(len(z),bool)
    for k,Z in active.items():
        ans &= ~Z[a*p**k+r//p**(L-k)]
    return ans
counts=[]
for dep in (5,4,3,2):
    L=6-dep
    Q=X[(X[:,:2*dep]==0).all(axis=1)&np.any(X[:,2*dep:2*dep+2],axis=1)]
    Z=np.zeros((p+1)*p**(L-1),dtype=bool)
    edge=np.unique(encode(Q[:,2*dep:]))
    edge=edge[allowed(edge,L)];Z[edge]=True
    while len(edge):
        nxt=[]
        for s in range(0,len(edge),200000):
            y=dec(edge[s:s+200000],L)
            for M in group:
                z=mm(y,M)%p
                assert not np.any(z[:,:2*dep]) and np.all(np.any(z[:,2*dep:2*dep+2],axis=1))
                k=encode(z[:,2*dep:])
                k=np.unique(k[(~Z[k])&allowed(k,L)])
                Z[k]=True;nxt.append(k)
        edge=np.concatenate(nxt)
    count=int(Z.sum());counts.append(count*p**dep);active[L]=Z
    print('depth',dep,'shortest prefix count',count,'covered profiles',counts[-1])
total=d['total'];mixed=d['mixed'];base=sum(d['counts'])+mixed
assert 26**46<16**46*p**6 and sum(counts)+mixed <= total
assert sum(counts)+mixed==4215171877
assert sum(counts)+mixed-base==1318136808
print('Certified exclusions',sum(counts)+mixed,'of',total,'increment',sum(counts)+mixed-base)
