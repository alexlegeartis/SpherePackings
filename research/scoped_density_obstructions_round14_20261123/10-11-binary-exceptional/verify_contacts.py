import itertools as I
import math
import numpy as np
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

def integer_shell(n, edges, hi):
    G = 2*np.eye(n, dtype="int64")
    for i,j in edges:
        G[i,j] = G[j,i] = -1
    Gs = Matrix(G.tolist())
    inv = Gs.inv()
    bounds = [math.isqrt(int(hi*inv[i,i])) for i in range(n)]
    seq = I.product(*(range(-u, u+1) for u in bounds))
    by = {}
    while True:
        v = np.asarray(list(I.islice(seq, 65536)), dtype="int64")
        if len(v)==0: break
        norms = np.sum((v@G)*v,axis=1)
        V=v[(norms>0)&(norms<=hi)]
        ids=(V%2)@np.array([1<<i for i in range(n)],dtype="int64")
        for x,s,d in zip(V,ids,norms[(norms>0)&(norms<=hi)]):
            by.setdefault((int(s),int(d)),[]).append(x)
    return Gs,by

for n in [4,6,7]:
    E = [(0,1),(1,2),(1,3)] if n==4 else [(i,i+1) for i in range(n-2)]+[(2,n-1)]
    G,b=integer_shell(n,E,6 if n==7 else 4)
    depth={0:0}
    for s,d in b:
        depth[s]=min(depth.get(s,d),d)
    assert len(depth)==2**n
    R = [Matrix(list(map(int,x))) for (s,d),xs in b.items() if d==2 for x in xs]
    assert len(R)=={4:24,6:72,7:126}[n]
    forms = Matrix([[int(r[i]*r[j]) for i in range(n) for j in range(i,n)] for r in R])
    assert forms.rank() == n*(n+1)//2
    assert sum((r*r.T for r in R),Matrix.zeros(n))==2*len(R)*G.inv()/n
    S={s for s,d in depth.items() if d==4}
    assert len(S)=={4:3,6:27,7:63}[n]
    A={s:Matrix(np.asarray(b[s,4],dtype="int64").T.tolist()) for s in S}
    assert {M.rank() for M in A.values()}==({4} if n==4 else {n-1})
    H={s:(M[:,1:]-M[:,0]*Matrix.ones(1,M.cols-1))/2 for s,M in A.items()}
    assert all(all(q.is_Integer for q in M) for M in H.values())
    H={s:Matrix([[int(x) for x in row] for row in M.tolist()]) for s,M in H.items()}
    if n==4:
        assert all(hermite_normal_form(M)==Matrix.eye(4) for M in H.values())
    else:
        lines={s:tuple(A[s].T.nullspace()[0]) for s in S}
        assert len(set(lines.values()))==len(lines)
        for a,b in I.combinations(S,2):
            if (a^b) in S:
                assert hermite_normal_form(Matrix.hstack(H[a],H[b]))==Matrix.eye(n)
                assert Matrix.hstack(*[Matrix(lines[s]) for s in (a,b,a^b)]).rank()==(2 if n==6 else 3)
    print(n,"OK",flush=True)
