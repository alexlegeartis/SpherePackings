"""Exact fifth-jet finite short-witness screen using the preceding packet."""
import gzip,json,runpy,sys
from pathlib import Path
v3=runpy.run_path(sys.argv.pop(1) if len(sys.argv)>1 else str(Path(__file__).resolve().parent.parent.parent/'scoped_density_round19_20261124'/'92-mw-cusp-index3'/'verify.py'))
g,np,p,mm,eq=(v3[x] for x in ('g','np','p','mm','eq'))
G,W,T,H,A,R,J,E,I=[v3[x] for x in ('G','W','T','H','A','R','J','E','I')]
S,n,D,invmod=[v3[x] for x in ('S','n','D','invmod')]
# Each subtraction is a residue section. Retain mod p^4 until division,
# so selected higher p-adic digits of H and row coordinates are preserved.
X=I.copy();jets=[]
for k in range(5):
    a=mm(X,A)%p;jets.append(a)
    X=mm((X-mm(a,R))%(p**(5-k)),(H+p**(5-k)//2)%(p**(5-k))-p**(5-k)//2)%(p**(5-k))
    assert not np.any(X%p)
    X=X//p
J3=np.column_stack(jets)
eq(J3[:,:6],v3['J3'])
eq(mm(E,J3)%p,np.column_stack((np.zeros_like(A),J3[:,:8])))
assert g['rank'](J3,p)==10
# In dimension 8 find inverse of independent rows; all quotient equations
# are checked on the 92 integral coordinate basis rows.
ids=[]
for i in range(n):
    if g['rank'](J3[ids+[i]],p)>len(ids):
        ids.append(i)
        if len(ids)==10:break
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
deep=np.array(D['depth4'],dtype='i8')
assert g['q'](deep)==24
assert not np.any(mm(deep,J3)[:8]%p)
assert np.any(mm(deep,J3)[8:]%p)
# All direct short-vector witnesses; no preliminary search needed.
seed0=[deep]
old=json.loads(Path(__file__).with_name('depth3.json').read_text())
old['depth']=3
literal=[old]
literal+=json.loads(gzip.decompress(Path(__file__).with_name('fifth_all.json.gz').read_bytes()))
for row in literal:
    w=np.array(row['w'],dtype='i8')
    z=mm(w,J3)%p
    d=row['depth']
    assert 0<=d<5 and not np.any(z[:2*d]) and np.any(z[2*d:2*d+2])
    assert g['q'](w)==row['norm'] and 16<=row['norm']<=24
    seed0.append(w)
X=mm(np.array(seed0,dtype='i8'),J3)%p
assert np.all(np.any(X,axis=1))
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
    assert not np.any(axis*ds[0])
    return axis*p**L+sum(ds[i]*p**(L-1-i) for i in range(L))
def dec(z,L):
    v=np.zeros((len(z),5,2),dtype='i8')
    axis=z//p**L; z=z%p**L;q=np.arange(len(z));k=5-L
    v[q,k,axis]=1
    for j in range(L):
        v[q,k+j,1-axis]=z//p**(L-1-j)%p
    return v.reshape(-1,10)
# Close shortest prefixes first; subtrees already covered are invariant
# under the normalizer and cannot contribute another exclusion.
active={}
def unblocked(z,L):
    axis=z//p**L;v=z%p**L
    ok=np.ones(len(z),bool)
    for k,covered in active.items():
        ok &= ~covered[axis*p**k+v//p**(L-k)]
    return ok
counts=[]
for depth in (4,3,2):
    L=5-depth
    Q=X[(np.count_nonzero(X[:,:2*depth],axis=1)==0)&np.any(X[:,2*depth:2*depth+2],axis=1)]
    Z=np.zeros((p+1)*p**(L-1),dtype=np.bool_)
    frontier=np.unique(encode(Q[:,2*depth:]))
    frontier=frontier[unblocked(frontier,L)]
    print('depth',depth,'initial residual prefixes',len(frontier),flush=True)
    Z[frontier]=True
    while len(frontier):
        rr=[]
        for i in range(0,len(frontier),200000):
            y=dec(frontier[i:i+200000],L)
            for m in group:
                vals=mm(y,m)%p
                assert not np.any(vals[:,:2*depth])
                z=encode(vals[:,2*depth:])
                z=z[(~Z[z])&unblocked(z,L)]
                z=np.unique(z).astype('i4')
                Z[z]=True
                if len(z):rr.append(z)
        frontier=np.concatenate(rr) if rr else []
    count=int(Z.sum())
    counts.append(count*p**depth)
    if depth==4:
        print('Leading directions:',np.flatnonzero(Z).tolist(),flush=True)
    print('depth',depth,'shortest witnessed prefixes',count,'excluded',counts[-1],flush=True)
    active[L]=Z
assert 24**46<16**46*47**5<26**46
mixed=(p+1)*(p**2+1)
total=(p+1)*p**4+mixed
print('Primitive excluded',sum(counts),'mixed excluded',mixed,'total',sum(counts)+mixed,'/',total,flush=True)
# The uncovered portion of this direct screen is a union of
# length-three prefixes, each leaving p^2 unrestricted lifts.
keys=np.arange((p+1)*p**2,dtype='i8')
leftover=keys[(~active[1][keys//p**2]) & (~active[2][keys//p]) & (~active[3][keys])]
assert p**2*len(leftover)==total-(sum(counts)+mixed)
pending=set(map(int,leftover))
orbits=[]
while pending:
    start=min(pending)
    seen={start}; stack=[start]
    while stack:
        y=dec(np.array([stack.pop()],dtype='i8'),3)
        for mat in group:
            z=int(encode((mm(y,mat)%p)[:,4:])[0])
            if z not in seen:
                seen.add(z);stack.append(z)
    assert seen<=pending
    pending-=seen
    orbits.append((len(seen), dec(np.array([start],dtype='i8'),3)[0].tolist()))
print('Uncovered primitive prefix orbits:',orbits)
