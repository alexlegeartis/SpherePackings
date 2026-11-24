"""Fourth-jet finite short-witness screen. Run in the index-47^3 packet folder."""
import runpy
from pathlib import Path
v3=runpy.run_path(str(Path(__file__).with_name('verify.py')))
g,np,p,mm,eq=(v3[x] for x in ('g','np','p','mm','eq'))
G,W,T,H,A,R,J,E,I=[v3[x] for x in ('G','W','T','H','A','R','J','E','I')]
S,n,D,invmod=[v3[x] for x in ('S','n','D','invmod')]
# Each subtraction is a residue section. Retain mod p^4 until division,
# so selected higher p-adic digits of H and row coordinates are preserved.
X=I.copy();jets=[]
for k in range(4):
    a=mm(X,A)%p;jets.append(a)
    X=mm((X-mm(a,R))%(p**(4-k)),H%(p**(4-k)))%(p**(4-k))
    assert not np.any(X%p)
    X=X//p
J3=np.column_stack(jets)
eq(J3[:,:6],v3['J3'])
eq(mm(E,J3)%p,np.column_stack((np.zeros_like(A),J3[:,:6])))
assert g['rank'](J3,p)==8
# In dimension 8 find inverse of independent rows; all quotient equations
# are checked on the 92 integral coordinate basis rows.
ids=[]
for i in range(n):
    if g['rank'](J3[ids+[i]],p)>len(ids):
        ids.append(i)
        if len(ids)==8:break
Z=invmod(J3[ids],p)
group=[]
for lab,k in [('U',5),('C',1),('F',1)]:
    X=g[lab]
    Tk=I.copy()
    for j in range(k):Tk=mm(Tk,T)
    eq(mm(X,T),mm(Tk,X))
    q=mm(Z,mm(X[ids],J3))%p
    eq(mm(X,J3)%p,mm(J3,q)%p)
    group.append(q)
seed0=np.vstack([v3['seed0'],v3['v'],v3['tail']]);d=v3['d']+[20,18]

seed=seed0.copy()
P=[mm(seed,J3)%p]
# Isometries prove image norms. The loop retains all 47 images, no heuristic sampling.
SX=mm(S,J3)%p
seed=seed%p
for k in range(p):
    P.append(mm(seed,SX)%p)
    seed=mm(seed,T)%p

# Build short literal sums/differences from the 121 checked input rows.
# Tested signed dot thresholds are exact; no claim of shell completeness.
nt=np.array(d,dtype='i8')
SG=mm(seed0,G)
pairs=[]
used=set()
sources=[I]
for i in range(3):
    sources.append(mm(sources[-1],g['U']))
sources.extend([g['C'],mm(g['C'],g['C']),g['F']])
Tpow=[I]
for k in range(p-1):Tpow.append(mm(Tpow[-1],T))
ops=[I]
for M in sources:
    ops.extend(mm(mm(M,Tt),S) for Tt in Tpow)
for O in ops:
    X=mm(seed0,O)
    dots=mm(X,SG.T)
    for sign in (1,-1):
        norms=nt[:,None]+nt[None,:]+sign*2*dots
        i,j=np.nonzero((norms>=16)&(norms<=22))
        for a,b in zip(i,j):
            w=X[a]+sign*seed0[b]
            code=tuple(w if w[next(i for i in range(n) if w[i])]>0 else -w)
            if code not in used:
                used.add(code)
                assert g['q'](w)==int(norms[a,b])
                pairs.append(w)
print('Short signed image pairs',len(pairs),flush=True)
for i0 in range(0,len(pairs),256):
    X=np.array(pairs[i0:i0+256],dtype='i8')
    P.append(mm(X,J3)%p)
    Xt=X%p
    for k in range(p):
        P.append(mm(Xt,SX)%p)
        Xt=mm(Xt,T)%p

X=np.unique(np.vstack(P),axis=0)
nz=np.any(X,axis=1)
assert len(pairs)==6643 and len(X)==317147
assert np.all(nz)
iv=np.array([0]+[pow(i,-1,p) for i in range(1,p)],dtype='i8')
def encode(y):
    L=y.shape[1]//2
    y=y.reshape(-1,L,2)
    axis=(y[:,0,0]==0).astype('i8')
    idx=np.arange(len(y)); a=y[idx,:,axis];b=y[idx,:,1-axis]
    i0=iv[a[:,0]]
    ds=[]
    for k in range(L):
        ds.append((i0*(b[:,k]-sum(ds[j]*a[:,k-j] for j in range(k))))%p)
    return axis*p**L+sum(ds[i]*p**(L-1-i) for i in range(L))
def dec(z,L):
    v=np.zeros((len(z),4,2),dtype='i8')
    axis=z//p**L; z=z%p**L;q=np.arange(len(z));k=4-L
    v[q,k,axis]=1
    for j in range(L):
        v[q,k+j,1-axis]=z//p**(L-1-j)%p
    return v.reshape(-1,8)
active=[]
for dept in range(4):
    Q=X[(np.count_nonzero(X[:,:2*dept],axis=1)==0)&np.any(X[:,2*dept:2*dept+2],axis=1)]
    L=4-dept
    Z=np.zeros(2*p**L,dtype=np.bool_)
    frontier=np.unique(encode(Q[:,2*dept:]))
    Z[frontier]=True
    while len(frontier):
        y=dec(frontier,L)
        rr=[]
        for m in group:
            vals=mm(y,m)%p
            assert not np.any(vals[:,:2*dept])
            z=encode(vals[:,2*dept:])
            zz=np.unique(z[~Z[z]])
            Z[zz]=True
            rr.append(zz)
        frontier=np.concatenate(rr)
    active.append(Z)
    print('depth',dept,'witness prefixes',np.count_nonzero(Z),flush=True)
ans=0
for d in (3,2,1,0):
    L=4-d
    Y=np.flatnonzero(active[d])
    axis=Y//p**L; v=Y%p**L
    ok=np.ones(len(Y),dtype=np.bool_)
    for e in range(d+1,4):
        k=4-e
        z=axis*p**k+v//p**(L-k)
        ok&=(~active[e][z])
    c=int(ok.sum())*p**d
    ans+=c
    print('depth',d,'incremental excluded',c,flush=True)
assert ans==3663078
print('Primitive index-p^4 profiles excluded',ans,'/',48*p**3,flush=True)


# Mixed (1,3) profiles are D times primitive index-p^2 lines.
# Shift depths by one; a depth-one vector supplies two digits and
# a depth-two vector supplies one digit.  No depth-three witness
# occurs among these seeds.
mixed=np.zeros(2*p*p,dtype=np.bool_)
Z=np.flatnonzero(active[1]); a=Z//p**3; b=Z%p**3
mixed[a*p*p+b//p]=True
Z=np.flatnonzero(active[2]); a=Z//p**2; b=Z%p**2
mixed[((a*p+b//p)[:,None]*p+np.arange(p)[None,:]).ravel()]=True
num=int(mixed.sum())
assert num==1968 and ans+num+1==3665047
print('Mixed (1,3) profiles excluded',num,'/',48*p)
# (2,2) is unique and contains the verified norm-20 row in D²L.
assert g['q'](v3['seed0'][-2])==20
assert 22**46<16**46*47**4
print('Total index-p^4 sublattice profiles excluded',ans+num+1,'/',48*p**3+48*p+1)
