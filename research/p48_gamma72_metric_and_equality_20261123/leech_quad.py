from itertools import combinations
from random import Random
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

def leechdata():
    # Standard extended binary Golay model (raw Leech norm multiplied by 8).
    g=sum(1<<j for j in (0,2,4,5,6,10,11))
    G=[0]
    for i in range(12):
        a=g<<i
        a^=(a.bit_count()%2)<<23
        G += [v^a for v in G]
    oc=[tuple(j for j in range(24) if z>>j&1) for z in G if z.bit_count()==8]
    assert len(set(G))==4096 and len(oc)==759
    return G,oc

def witness():
    G,oc=leechdata()
    seed=Random(20261123)
    def row():
        t=seed.randrange(3)
        if t==0:
            i,j=seed.sample(range(24),2)
            v=[0]*24
            v[i],v[j]=seed.choice((4,-4)),seed.choice((4,-4))
        elif t==1:
            v=[0]*24
            arr=seed.choice(oc)
            m=seed.randrange(128);m|=(m.bit_count()%2)<<7
            for i,j in enumerate(arr):v[j]=2-4*((m>>i)&1)
        else:
            m=seed.choice(G);j=seed.randrange(24)
            v=[1-2*((m>>i)&1) for i in range(24)]
            v[j]*=-3
        assert sum(x*x for x in v)==32
        return v
    rows=[row() for _ in range(80)]
    B=hermite_normal_form(Matrix(rows).T).T
    assert B.shape==(24,24) and abs(B.det())==2**36
    U=8**12*B.inv()
    assert all(x.q==1 for x in U)
    # All evaluation rows are on this lattice's integer-coordinate mod-2 classes.
    Gram=B*B.T/8
    assert all(x.q==1 for x in Gram) and all(Gram[i,i]%2==0 for i in range(24))
    assert Gram.det()==1
    pairs=[(i,j) for i in range(24) for j in range(i,24)]
    Q=sum((int((Gram[i,i]/2 if i==j else Gram[i,j]))%2)<<k for k,(i,j) in enumerate(pairs))
    assert Q!=0
    UU=[[int(x) for x in u] for u in U.tolist()]
    pivots={}
    for attempt in range(1,3001):
        raw=row()
        zz=[sum(x*UU[i][j] for i,x in enumerate(raw)) for j in range(24)]
        assert all(x%2**36==0 for x in zz)
        z=[(x//2**36)&1 for x in zz]
        H=sum((z[i]&z[j])<<k for k,(i,j) in enumerate(pairs))
        assert (H&Q).bit_count()%2==0
        while H:
            k=H.bit_length()
            if k in pivots:H^=pivots[k]
            else:
                pivots[k]=H
                break
        if len(pivots)==299:
            print('PASS: Leech minimum-vector classes admit exactly one independent binary quadratic identity; rank',len(pivots),'after',attempt,'rows')
            return
    raise AssertionError(len(pivots))
if __name__=='__main__':witness()
