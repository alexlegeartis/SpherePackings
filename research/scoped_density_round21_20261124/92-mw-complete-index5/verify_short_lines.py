"""Exact short F_p-line orbits in the depth-two quotient B=K2/K5."""
import runpy
from fractions import Fraction
from pathlib import Path
here=Path(__file__).resolve().parent
closure=runpy.run_path(str(here/'verify_fifth_complete.py'))
base=closure['d']
p,np,mm,group=(base[k]for k in ('p','np','mm','group'))
X=np.vstack((base['X'],[row['jet']for row in closure['rows']]))
assert np.all(X[:,:4]==0)
iv=np.array([0]+[pow(i,-1,p) for i in range(1,p)],dtype='i8')

def encode(Y,L):
    y=Y[:,-2*L:].reshape(-1,L,2)
    axis=(y[:,0,0]==0).astype('i8')
    s=np.arange(len(y))
    scale=iv[y[s,0,axis]]
    assert np.all(scale)
    v=(y*scale[:,None,None])%p
    assert np.all(v[s,0,axis]==1) and not np.any(axis*v[s,0,1-axis])
    if L>=2:
        j=v[s,1,axis].copy()
        a=v[:,1].copy()
        if L==3:
            v[:,2]=(v[:,2]-j[:,None]*a+(j*(j-1)//2)[:,None]*v[:,0])%p
        v[:,1]=(a-j[:,None]*v[:,0])%p
        assert not np.any(v[s,1,axis])
    z=axis*p+v[s,0,1-axis]
    if L>=2:z=z*p+v[s,1,1-axis]
    if L==3:z=(z*p+v[s,2,axis])*p+v[s,2,1-axis]
    return z

def dec(z,L):
    m=1 if L==1 else (2 if L==2 else 4)
    axis=z//p**m
    q=z%p**m; rr=np.arange(len(z));d=5-L
    y=np.zeros((len(z),5,2),dtype='i8')
    y[rr,d,axis]=1; y[rr,d,1-axis]=q//p**(m-1)
    if L>=2:y[rr,d+1,1-axis]=q//p**(m-2)%p
    if L==3:y[rr,d+2,axis]=q//p%p; y[rr,d+2,1-axis]=q%p
    return y.reshape(-1,10)

counts=[]
for d in (2,3,4):
    L=5-d; m=1 if L==1 else (2 if L==2 else 4)
    Z=np.zeros((p+1)*p**(m-1),dtype=bool)
    Y=X[(X[:,:2*d]==0).all(axis=1)&np.any(X[:,2*d:2*d+2],axis=1)]
    queue=np.unique(encode(Y,L)); Z[queue]=True
    while len(queue):
        y=dec(queue,L); assert np.array_equal(encode(y,L),queue)
        next_=[]
        for M in group:
            w=mm(y,M)%p
            assert not np.any(w[:,:2*d]) and np.all(np.any(w[:,2*d:2*d+2],axis=1))
            z=encode(w,L)
            z=np.unique(z[~Z[z]]); Z[z]=True;next_.append(z)
        queue=np.concatenate(next_)
    counts.append(int(Z.sum()))
assert counts==[68724,828,12]
print('Canonical translation orbit counts by depth 2,3,4:',counts)
short=p*(counts[0]+counts[1])+counts[2]
assert short==3268956
assert 0<short<=(p**6-1)//(p-1)
print('Explicit short F_p-lines in K2/K5:',short)
def gauss(n,k):
    r=Fraction(1,1)
    for j in range(k):r*=Fraction(p**(n-j)-1,p**(k-j)-1)
    assert r.denominator==1
    return r.numerator
inc=short*gauss(5,1)
sec=inc+short*(short-1)
two=max(inc-short*(short-1)//2,(inc**2+sec-1)//sec)
assert two==10955198305026
assert 0<two<=gauss(6,2)
print('Two-dimensional subspaces containing a witnessed line: at least',two,'of',gauss(6,2))
C={3:(p+1)*p**2,4:(p+1)*p**2*(p**2+p)+(p+1)*p,5:(p**6-1)//(p-1)}
expected={5:9179147596354897233629066132555953562124,6:4184518740099439915762768085636498388,7:857658426082247424333241185192}
for d in (5,4,3):
    f=lambda e,k:k*gauss(4,d-e)*p**((6-e)*(d-e))
    t1=f(1,short);t2=f(2,two);t0=sum(f(e,C[e])for e in range(3,d+1))
    ans=t1+t2+t0
    assert ans==expected[10-d]
    assert ans<=gauss(10,d) and 24**46<16**46*p**(10-d)
    print('Index',f'{p}^{10-d}','short-witness sublattices at least',ans,'disjoint contributions',[t1,t2,t0])
