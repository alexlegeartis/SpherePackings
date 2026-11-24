"""Exact MW92 cusp duality replay, including one- and three-jet forms."""
from fractions import Fraction
from math import prod
from pathlib import Path
import gzip,json,sys
import numpy as np

here=Path(__file__).resolve().parent
f=Path(sys.argv[1]) if len(sys.argv)>1 else here.parent.parent/'scoped_density_round17_20261124'/'92-MW-cusp-jets'/'input.json.gz'
p=47
src=json.loads(gzip.decompress(f.read_bytes()))
T,G,A,U=[np.asarray(src[k],dtype='i8') for k in ('T','G','A','U')]
n=len(G)
I=np.eye(n,dtype='i8')
D=I-T
def eq(A,B):
    assert np.array_equal(A,B)
def mul(X,Y):
    assert X.shape[-1]==Y.shape[0]
    assert X.shape[-1]*int(np.max(np.abs(X),initial=0))*int(np.max(np.abs(Y),initial=0))<2**62
    return X@Y
def powm(X,k):
    R=I.copy()
    for _ in range(k):R=mul(R,X)
    return R
def pow_mod(X,k,m):
    Y=I.copy()
    for _ in range(k):Y=mul(Y,X)%m
    return Y
def det_bareiss(g):
    a=[[int(v) for v in row] for row in g]
    old=1
    for i in range(len(a)-1):
        pivot=a[i][i]; assert pivot>0
        for j in range(i+1,len(a)):
            for k in range(i+1,len(a)):
                q,r=divmod(a[j][k]*pivot-a[j][i]*a[i][k],old)
                assert not r
                a[j][k]=q
        old=pivot
    return a[-1][-1]
eq(G,G.T);assert det_bareiss(G)==p**14
Z=np.zeros_like(I);S=Z.copy();H=Z.copy();B=Z.copy();Q=Z.copy();TT=I.copy()
for j in range(p):
    S+=TT;H-=j*TT
    if pow(j,(p-1)//2,p)==1:B+=TT
    elif j:B-=TT
    if j%2:Q-=TT  # inverse of 1+T
    TT=mul(TT,T)
eq(S,Z);eq(TT,I);eq(mul(D,H),p*I)
eq(mul(mul(T,G),T.T),G)
eq(mul(I+T,Q),I)
D7=powm(D,7)
def solve(M,X,m):
    A=np.column_stack([M%m,X%m]).astype('i8')
    assert (m-1)**2+m<2**62
    for i in range(n):
        j=next(j for j in range(i,n) if A[j,i]%m)
        A[[i,j]]=A[[j,i]]
        A[i]=A[i]*pow(int(A[i,i]),-1,m)%m
        col=A[:,i].copy();col[i]=0
        A=(A-col[:,None]*A[i])%m
    eq(A[:,:n],I)
    R=A[:,n:]
    return np.where(R>m//2,R-m,R)
Y,V=np.hsplit(solve(G,np.column_stack((D7,p*I)),16777213),2)
eq(mul(G,Y),D7);eq(mul(G,V),p*I)
assert (np.max(abs(D7)),np.max(abs(Y)),np.max(abs(V)))==(8556,80066,456438)
assert not np.any(V.diagonal()%2)
print('Dual certificate: D^7 L* = L, and even Gram 47 G^{-1} for D^16 L.')
J=powm(U,23)
eq(mul(J,J),-I)
eq(mul(mul(J,G),J.T),G)
eq(mul(J,T),mul(powm(T,46),J))
eq(mul(J,B),-mul(B,J))
eq(mul(B,B),-p*I);eq(mul(mul(B,G),B.T),p*G)
pi=mul(D,Q)
eq(mul(pi,J),-mul(J,pi))
eq(mul(pi,G),-mul(G,pi.T))
K=mul(J,B)
eq(mul(K,K),-p*I)
eq(mul(mul(K,G),K.T),p*G)
eq(mul(K,G),-mul(G,K.T))
eq(mul(J,A)%p,mul(A,np.asarray([[0,1],[-1,0]],'i8'))%p)
print('J=U^23 reverses T, J²=-I; skew uniformizer and the common 47-modular similarity verified.')

eq(mul(pow_mod(D,38,p),G)%p,(4*mul(A,A.T))%p)
i,j=next((i,j)for i in range(n)for j in range(i) if (A[i,0]*A[j,1]-A[i,1]*A[j,0])%p)
a,b,c,d=map(int,A[[i,j]].ravel())
R=np.array([[d,-b],[-c,a]],dtype='i8')*pow((a*d-b*c)%p,-1,p)@I[[i,j]]%p
# Build the jets on each integral basis row independently.
Z0=I.copy();js=[]
for k in range(3):
    a=mul(Z0,A)%p;js.append(a)
    t=p**(3-k)
    Z0=mul((Z0-mul(a,R))%t,H%t)%t
    assert not np.any(Z0%p);Z0//=p
J3=np.column_stack(js)
F=np.array([[43,16,23,32,4,0],[39,45,15,23,0,4],[20,15,43,0,0,0],
            [32,20,0,43,0,0],[4,0,0,0,0,0],[0,4,0,0,0,0]],dtype='i8')
eq(mul(pow_mod(D,36,p),G)%p,mul(mul(J3,F),J3.T)%p)
R3=np.vstack([R%p,mul(R,D)%p,mul(mul(R,D)%p,D)%p])
eq(mul(R3,J3)%p,np.eye(6,dtype='i8'))
Jbar=mul(mul(R3,J),J3)%p
Ubar=mul(mul(R3,U),J3)%p
def orientation(x):
    x=x.reshape(3,2)
    a,b=x.T
    assert a[0]
    r=[]
    for i in range(3):r.append(int((b[i]-sum(r[j]*a[i-j] for j in range(i)))*pow(int(a[0]),-1,p)%p))
    return tuple(r)
for u,v in [((2,27,21),(23,15,14)),((6,24,0),(39,18,29))]:
    x=np.asarray([1,u[0],0,u[1],0,u[2]],'i8')
    y=np.asarray([1,v[0],0,v[1],0,v[2]],'i8')
    shift=lambda z:np.r_[0,0,z[:4]]
    assert all(int(y@F@z)%p==0 for z in [x,shift(x),shift(shift(x))])
    assert orientation(x@Jbar%p)==orientation(x@Ubar%p)==v
    print('Three-jet pairing',u,'<->',v)
print('Exact duality identities and unresolved prefix partners verified.')
