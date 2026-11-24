"""Exact replay for translated cusp-jet obstructions; NumPy uses guarded integers."""
import json,runpy,sys
from pathlib import Path
here=Path(__file__).resolve().parent
base=Path(sys.argv[1]) if len(sys.argv)>1 else here.parent.parent/'scoped_density_round17_20261124'/'92-MW-cusp-jets'/'verify.py'
r18=Path(sys.argv[2]) if len(sys.argv)>2 else here.parent.parent/'scoped_density_round18_20261124'/'92-mw-cusp-index2'
g=runpy.run_path(str(base))
np,p,mm,eq=g['np'],g['p'],g['mm'],g['eq']
G,W,T,H,A,R,J,E,I=[g[x] for x in ('G','W','T','H','A','R','J','E','I')]
D,ht,sections,trim,sub,mul,power=g['D'],g['ht'],g['sections'],g['trim'],g['sub'],g['mul'],g['power']
n=len(G)
# x(T) -> T^16 x(-1/T), y likewise at degree 24.
def reverse(x,deg):
    return trim([g['amin'](0,x[deg-i]) if (deg-i)%2 else x[deg-i] for i in range(deg+1)])
Q=[]
for x,y in sections:
    X,Y=reverse(x,16),reverse(y,24)
    assert sub(power(Y,2),power(X,3))==g['f']
    Q.append((X,Y))
B=np.array([[ht(x,y) for y in sections] for x in Q],dtype='i8')
def invmod(A,m):
    A=A.copy()%m
    d=len(A);B=np.eye(d,dtype='i8')
    for j in range(d):
        i=next(i for i in range(j,d) if A[i,j])
        A[[i,j]]=A[[j,i]];B[[i,j]]=B[[j,i]]
        k=pow(int(A[j,j]),-1,m);A[j]=A[j]*k%m;B[j]=B[j]*k%m
        t=A[:,j].copy();t[j]=0
        A=(A-t[:,None]*A[j])%m;B=(B-t[:,None]*B[j])%m
    eq(A,np.eye(d,dtype='i8'))
    return B
m=65521
S=mm(mm(invmod(W,m),B)%m, invmod(mm(G,W.T),m))%m
S=np.where(S>m//2,S-m,S)
eq(mm(mm(mm(W,S),G),W.T),B)
eq(mm(mm(S,G),S.T),G);eq(mm(S,S),I)
print('Cusp inversion matrix reconstructed and checked exactly; max abs entry',np.max(abs(S)))
# Each subtraction is a residue section. Retain mod p^3 until division,
# so selected higher p-adic digits of H and row coordinates are preserved.
X=I.copy();jets=[]
for k in range(3):
    a=mm(X,A)%p;jets.append(a)
    X=mm((X-mm(a,R))%(p**(3-k)),H%(p**(3-k)))%(p**(3-k))
    assert not np.any(X%p)
    X=X//p
J3=np.column_stack(jets)
eq(J3[:,:4],J)
eq(mm(E,J3)%p,np.column_stack((np.zeros_like(A),J)))
assert g['rank'](J3,p)==6
# In dimension 6 find inverse of independent rows; all quotient equations
# are checked on the 92 integral coordinate basis rows.
ids=[]
for i in range(n):
    if g['rank'](J3[ids+[i]],p)>len(ids):
        ids.append(i)
        if len(ids)==6:break
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
def pref(x):
    x=np.array(x).reshape(3,2).tolist()
    d=next((i for i,w in enumerate(x) if any(w)),3)
    if d==3:return None
    rows=x[d:]
    axis=int(not rows[0][0]);f=pow(rows[0][axis],-1,p)
    ans=[]
    for r,row in enumerate(rows):
        ans.append(f*(row[1-axis]-sum(ans[j]*rows[r-j][axis] for j in range(r)))%p)
    return (axis,*ans)
def row(c):
    axis,*ans=c
    x=np.zeros((3,2),dtype='i8');d=3-len(ans)
    x[d,axis]=1
    x[d:,1-axis]=ans
    return x.ravel()
def prefixes(x):
    seen=set(y for z in x if (y:=pref(z)) is not None)
    stack=list(seen)
    while stack:
        c=stack.pop();v=row(c)
        for Q in group:
            z=pref(mm(v,Q)%p); assert z is not None and len(z)==len(c)
            if z not in seen:
                seen.add(z);stack.append(z)
    return seen
w18=json.loads((r18/'witnesses.json').read_text())
w20=json.loads((r18/'norm20.json').read_text())
seed=np.array(D['q16rows']+D['q18rows']+[v[-1] for v in w18]+D['W']+[D['depth2'],w20],dtype='i8')
d=[16]*17+[18]*2+[v[-2] for v in w18]+[16]*92+[20]*2
assert len(seed)==119 and [g['q'](x) for x in seed]==d
assert np.all(mm(seed[-2:],J3)[:,:4]%p==0)
seed0=seed.copy()
P=[mm(seed,J3)%p]
# Isometries prove image norms. The loop retains all 47 images, no heuristic sampling.
SX=mm(S,J3)%p
seed=seed%p
for k in range(p):
    P.append(mm(seed,SX)%p)
    seed=mm(seed,T)%p

# Build short literal sums/differences from the same 119 pinned rows.
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
        i,j=np.nonzero((norms>=16)&(norms<=20))
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
Z=prefixes(np.vstack(P))
print('Number of witnessed prefixes of depth 2,1,0:',*[sum(len(v)==i for v in Z) for i in (2,3,4)])
remaining=[]
counts=[0,0,0]
for a in (0,1):
    for x in (range(p) if a==0 else [0]):
        for y in range(p):
            for z in range(p):
                key=(a,x,y,z)
                for k in (2,3,4):
                    if key[:k] in Z:
                        counts[k-2]+=1
                        break
                else:
                    remaining.append(key)
print('Covered at depth 2,1,0 without duplicates:',counts)
print('Exclusions including 48 mixed elementary divisor profiles:',sum(counts)+48)
print('Remaining profiles:',len(remaining),remaining)
assert len(pairs)==1460
assert counts==[26508,12972,66528] and len(remaining)==24
assert 20**46<16**46*47**3

# A single literal norm-18 seed found by cusp inversion of a previous search
# orbit. The large discovery archive is not a replay input.
tail=np.array([-16,-32,-3,22,17,-9,11,-23,-8,24,21,12,22,-28,30,18,-20,19,-2,20,20,-2,8,-3,22,13,-6,14,8,-8,6,0,-14,-1,-6,4,9,6,6,-10,-9,1,-10,6,-5,-6,-4,-5,-2,9,2,11,-1,-6,-1,-6,-2,0,4,0,1,-1,1,-4,6,6,2,-4,-1,2,2,1,1,2,1,4,1,-1,0,0,2,2,0,1,1,-1,0,-1,-2,1,-1,-1],dtype='i8')
assert g['q'](tail)==18
tq=mm(tail,J3)%p
assert pref(tq)==(0,20,12,37)
new=prefixes(np.array([tq]))
assert len(new)==12
remaining=[x for x in remaining if x not in new]
print('After literal final-orbit witness:',106080-len(remaining),'remaining',remaining)
assert len(remaining)==12
assert prefixes(np.array([row(remaining[0])]))==set(remaining)

# Final short witness, found by matching three-jet syndromes and taking
# a difference of two known norm-18 rows. Check the resulting integer row
# directly; discovery inputs are not needed.
v=np.array([3,24,-16,-24,-19,2,-3,-11,-10,-10,-2,-13,-21,37,-6,0,18,-13,-5,-23,-19,8,3,2,-15,-9,10,-23,-23,2,9,16,-5,-3,4,9,-6,6,2,3,3,2,-7,-2,10,-6,-5,-4,0,-3,-4,2,-4,2,2,1,3,1,-1,1,1,-2,4,1,-1,0,1,-1,1,0,-3,-2,2,1,1,-2,0,-1,1,0,-1,3,-2,-1,1,1,1,0,1,-1,-1,-2],dtype='i8')
assert g['q'](v)==20
vq=mm(v,J3)%p
assert pref(vq)==(0,3,12)
new=prefixes(np.array([vq]))
assert all(x[:3] in new for x in remaining)
print('All 106080 translation-invariant index-47^3 sublattices contain a witness of squared norm <=20.')
