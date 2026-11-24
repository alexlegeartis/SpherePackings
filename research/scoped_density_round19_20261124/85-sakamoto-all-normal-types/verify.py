"""Exact replay for one chart of every admissible low-determinant normal Gram type, for each Sakamoto parent."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations_with_replacement
from pathlib import Path
import gzip, hashlib, json, sys
import numpy as np


def maxabs(x):
    return max(int(x.max(initial=0)), -int(x.min(initial=0)))


def mm(a,b):
    assert a.shape[-1] == b.shape[0] and a.shape[-1]*maxabs(a)*maxabs(b)<2**62
    return a@b


def dotrows(a,b):
    assert a.shape==b.shape and a.shape[-1]*maxabs(a)*maxabs(b)<2**62
    return np.einsum('ij,ij->i',a,b)


def gram_norm(a,G):
    return dotrows(mm(a,G),a)


def proposal_hits(a,b,target):
    """Fast exact BLAS dot check: all products and all possible partial sums fit in 53 bits."""
    assert a.shape[1]==b.shape[1] and a.shape[1]*maxabs(a)*maxabs(b)<2**52
    return np.where(a.astype(float)@b.T.astype(float)==target)


def expand(a,M5,M23):
    out=[]
    for i in range(5):
        b=a
        for j in range(23):
            out.append(b)
            b=mm(b,M23)
        assert np.array_equal(a,b)
        a=mm(a,M5)
    assert max(maxabs(x) for x in out)<2**61
    return np.unique(np.concatenate(out+[-x for x in out]),axis=0)


def stock(G,M5,M23,A):
    # Discovery only: all retained vectors are independently certified below.
    rng=np.random.default_rng(365503)
    B=expand(A,M5,M23)
    assert len(B)==10120 and np.all(gram_norm(B,G)==8)
    old={bytes(b) for b in B}
    for rep in range(12):
        C=B[:88] if not rep else B[rng.choice(len(B),64,replace=False)]
        js,ks=proposal_hits(mm(C,G),B,4)
        assert max(maxabs(C),maxabs(B))<2**61
        W=np.unique(C[js]-B[ks],axis=0)
        hit=[w for w in W if bytes(w) not in old]
        np.random.default_rng(rep).shuffle(hit)
        seeds=[]
        for w in hit:
            if bytes(w) in old:
                continue
            V=expand(w[None,:],M5,M23)
            assert np.all(gram_norm(V,G)==8)
            old.update(bytes(v) for v in V)
            seeds.append(V)
            if len(seeds)>=125:
                break
        if not seeds:
            break
        B=np.unique(np.concatenate([B,*seeds]),axis=0)
        if len(B)>=200000:
            break
    assert np.all(gram_norm(B,G)==8)
    return B


def omega(A):
    best=0
    def col(P):
        order=[]
        bounds=[]
        d=0
        while P:
            d+=1
            Q=P
            while Q:
                b=Q&-Q
                u=b.bit_length()-1
                P^=b;Q^=b;Q&=~A[u]
                order.append(u);bounds.append(d)
        return order,bounds
    def rec(P,s):
        nonlocal best
        if not P:
            best=max(best,s);return
        order,ub=col(P)
        for u,k in zip(reversed(order),reversed(ub)):
            if s+k<=best:
                return
            rec(P&A[u],s+1)
            P^=1<<u
    rec((1<<len(A))-1,0)
    return best


def adjugate(g):
    p=g.tolist()
    def co(i,j):
        r=[u for u in range(3) if u!=j]
        c=[u for u in range(3) if u!=i]
        return (-1)**(i+j)*(p[r[0]][c[0]]*p[r[1]][c[1]]-p[r[0]][c[1]]*p[r[1]][c[0]])
    adj=np.array([[co(i,j) for j in range(3)] for i in range(3)],dtype='i8')
    D=int(mm(g,adj)[0,0])
    assert D>0 and np.array_equal(mm(g,adj),D*np.eye(3,dtype='i8'))
    return adj,D


def normalized(a,b,c):
    v=sorted(map(abs,(a,b,c)))
    if a*b*c<0:
        v[2]=-v[2]
    return tuple(v)


def types():
    out=set()
    for a,b,c in combinations_with_replacement(range(5),3):
        for z in ([c,-c] if a else [c]):
            if min(a+b+z,a-b-z,-a+b-z,-a-b+z)<-8:
                continue
            d=512+2*a*b*z-8*(a*a+b*b+z*z)
            if 0<d<=384:
                out.add((a,b,z))
    return out


def certify(g,BG,ns,C):
    adj,D=adjugate(g)
    E=mm(BG,C.T)
    F0=mm(E,adj)
    weights=np.array([1,D,D*D],dtype='i8')
    labels=mm(F0%D,weights)
    score=D*(8-ns)+dotrows(E,F0)
    # The adjugate-coordinate group is exactly image(adj mod D).
    V={tuple([0,0,0])}
    todo=list(V)
    while todo:
        u=todo.pop()
        for v in adj.tolist():
            q=tuple((a+b)%D for a,b in zip(u,v))
            if q not in V:
                todo.append(q);V.add(q)
    assert len(V)==D
    P=np.array(sorted(V),dtype='i8')
    Pkey=mm(P,weights)
    B={}
    # Group maxima by explicit array reduction, then check every winner's literal certificate.
    o=np.argsort(labels,kind='stable')
    k=labels[o]
    starts=np.r_[0,1+np.flatnonzero(k[1:]!=k[:-1])]
    maxima=np.maximum.reduceat(score[o],starts)
    assert set(Pkey)>=set(k[starts])
    for residue,h in zip(k[starts],maxima):
        if h<=0 or residue==0:
            continue
        # Negation is equally witnessed; keep this explicit.
        v=tuple(np.unravel_index(int(residue),(D,D,D))[::-1])
        u=sum(((-x)%D)*int(w) for x,w in zip(v,weights))
        B[int(residue)]=max(B.get(int(residue),0),int(h))
        B[u]=max(B.get(u,0),int(h))
    # Norm 8 in horizontal period is certified by a literal source norm-8 zero *normal*.
    # Do not confuse zero residue with zero normal label.
    assert np.any((ns==8)&np.all(E==0,axis=1))
    hcast=np.array([B.get(int(k),0) for k in Pkey],dtype='i8')
    dkeys=((P[:,None,:]-P[None,:,:])%D)@weights
    hpair=np.vectorize(lambda x:B.get(int(x),0),otypes=['i8'])(dkeys)
    delta2=F(0)
    def ratio2(M,h):
        return F(4*M*M*(8*D-h)**85,D*(8*D)**85)
    for h in sorted({0,*B.values()}):
        if h>=8*D:
            break
        if ratio2(D,h)<=F(16,25):
            break
        ids=np.flatnonzero((hcast<=h)&(Pkey!=0))
        b=hpair[ids[:,None],ids[None,:]]<=h
        A=[sum(1<<int(j) for j in np.flatnonzero(row) if i!=j) for i,row in enumerate(b)]
        m=1+omega(A)
        assert ratio2(m,h)<=F(16,25),(D,h,m,ratio2(m,h))
        delta2=max(delta2,ratio2(m,h))
    return D,len(B),delta2


def main():
    path=Path(__file__).with_name('input.json.gz')
    raw=json.loads(gzip.decompress(path.read_bytes()))
    sha={'L':'f95d33297e2713f696b13d5e123d24b86f2dd3212ee514fea6cafc2548f11c06',
         'L0':'53cf44a2d5ab1af9b4c2d93ad8b54d84041f1a4871874e6d8be4ffcce2bb247a'}
    allowed=types()
    assert len(allowed)==26
    for tag in sys.argv[1:] or ['L0','L']:
        d=raw[tag]
        assert hashlib.sha256(json.dumps(d['G'],separators=(',',':')).encode()).hexdigest()==sha[tag]
        G,M5,M23=[np.array(d[q],dtype='i8') for q in ('G','M5','M23')]
        if not np.array_equal(mm(mm(M5,G),M5.T),G):M5=M5.T
        if not np.array_equal(mm(mm(M23,G),M23.T),G):M23=M23.T
        assert np.array_equal(mm(mm(M5,G),M5.T),G)
        assert np.array_equal(mm(mm(M23,G),M23.T),G)
        B=stock(G,M5,M23,np.array(d['seeds'],dtype='i8'))
        rng=np.random.default_rng(1200+len(tag))
        a=B[rng.choice(len(B),180,replace=False)]
        j,k=proposal_hits(mm(a,G),B,3)
        N=np.unique(a[j]-B[k],axis=0)
        assert np.all(gram_norm(N,G)==10)
        BG=mm(np.concatenate([B,N]),G)
        ns=np.r_[np.full(len(B),8,dtype='i8'),np.full(len(N),10,dtype='i8')]
        del a,N,B
        cs=set()
        for t,C in d['charts']:
            C=np.array(C,dtype='i8')
            g=mm(mm(C,G),C.T)
            assert tuple(g.diagonal())==(8,8,8)
            shape=normalized(int(g[0,1]),int(g[0,2]),int(g[1,2]))
            assert tuple(t)==shape and shape not in cs
            cs.add(shape)
            D,no,delta2=certify(g,BG,ns,C)
            assert 0<D<=384
            print(tag,shape,'det',D,'blocked',no,'exact cap squared ≤16/25',float(delta2),flush=True)
        assert cs==allowed
        print('PASS',tag,len(cs),flush=True)
if __name__=='__main__':
    main()
