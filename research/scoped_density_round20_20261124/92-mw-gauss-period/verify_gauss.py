"""Exact period-23 Gauss similarity and deep jet upper witnesses."""
import json,gzip,sys
from pathlib import Path
import numpy as np
p=47
here=Path(__file__).resolve().parent
pinned=Path(sys.argv[1]) if len(sys.argv)>1 else here.parent.parent/'scoped_density_round17_20261124'/'92-MW-cusp-jets'/'input.json.gz'
obj=json.loads(gzip.decompress(pinned.read_bytes()))
G,T=np.array(obj['G'],'i8'),np.array(obj['T'],'i8')
n=len(G)
I=np.eye(n,dtype='i8')
def mul(X,Y):
    assert X.shape[-1]==Y.shape[0]
    assert X.shape[-1]*int(np.max(np.abs(X),initial=0))*int(np.max(np.abs(Y),initial=0))<2**62
    return X@Y
H,TT,Gaus,Tot=[np.zeros_like(G) for _ in range(4)]
TT+=I
Leg=[0]+[1 if pow(j,(p-1)//2,p)==1 else -1 for j in range(1,p)]
for j in range(p):
    H-=j*TT
    Gaus+=Leg[j]*TT
    Tot+=TT
    TT=mul(TT,T)
D=I-T
assert np.array_equal(TT,I) and not np.any(Tot)
assert np.array_equal(mul(D,H),p*I)
assert np.array_equal(mul(mul(T,G),T.T),G)
assert np.array_equal(mul(Gaus,Gaus),-p*I)
assert np.array_equal(mul(mul(Gaus,G),Gaus.T),p*G)
for i in range(p):
    assert sum(Leg[j]*Leg[(j+i)%p] for j in range(p))==(p-1 if i==0 else -1)
print('Quadratic Gauss operator verified: B²=-47I, BGBᵀ=47G')
R=json.loads((here/'deep_upper_witnesses.json').read_text())
for d,q,v in R:
    v=np.asarray(v,'i8')
    assert v.shape==(n,) and np.any(v)
    w=v.copy()
    assert int(mul(mul(w,G),w))==q
    for i in range(d):
        w=mul(w,H);assert not np.any(w%p);w//=p
    for i in range(d):w=mul(w,D)
    assert np.array_equal(w,v)
    k=next(k for k in range(100) if q**46<16**46*p**k)
    print(f'Depth {d}: literal squared norm {q}; blocks every profile containing this depth at index {p}^{k} or larger')
