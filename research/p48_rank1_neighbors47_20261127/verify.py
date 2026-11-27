"""Exact harmonic moments, rational shell polytopes, and phase graphs."""
from fractions import Fraction as Q
from itertools import combinations
from math import gcd, isqrt, lcm
import sys

assert __debug__


def mul(a,b,n=3):
    return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(n)]


def power(a,k,n=3):
    z=[1]+[0]*(n-1)
    for _ in range(k):
        z=mul(z,a,n)
    return z


def ordinary_counts():
    n=6
    E=[1]+[240*sum(d**3 for d in range(1,k+1) if k%d==0)
           for k in range(1,n)]
    D=[0,1]+[0]*(n-2)
    for k in range(1,n):
        v=[1]+[0]*(n-1)
        v[k]=-1
        D=mul(D,power(v,24,n),n)
    u=power(E,6,n)
    v=mul(power(E,3,n),D,n)
    z=power(D,2,n)
    return [a-1440*b+125280*c for a,b,c in zip(u,v,z)]


def rref(T,n):
    T=[list(map(Q,r)) for r in T]
    piv=[]
    for j in range(n):
        k=len(piv)
        ix=next((i for i in range(k,len(T)) if T[i][j]),None)
        if ix is None:
            continue
        T[k],T[ix]=T[ix],T[k]
        a=T[k][j]
        T[k]=[x/a for x in T[k]]
        for i,r in enumerate(T):
            if i!=k and r[j]:
                a=r[j]
                T[i]=[x-a*y for x,y in zip(r,T[k])]
        piv.append(j)
    assert all(not any(r[n:]) for r in T[len(piv):])
    return piv,T[:len(piv)]


def Fourier_relations(k):
    basis=[]
    for c in range(k//12+1):
        for b in (0,1):
            a,remainder=divmod(k-12*c-6*b,4)
            if a<0 or remainder:
                continue
            z=mul([1,-72,2484],power([0,1,-24],c))
            z=mul(z,power([1,240,2160],a))
            z=mul(z,power([1,-504,-16632],b))
            basis.append(z+[0])
    piv,R=rref(basis,3)
    for j in range(3):
        if j not in piv:
            a=[Q(i==j) for i in range(3)]
            for i,p in enumerate(piv):
                a[p]=-R[i][j]
            yield a


def zonal(a,w,m,ell):
    coeff,z=Q(1),Q(0)
    for j in range(ell//2+1):
        k=ell-2*j
        if j:
            coeff*=Q(-(k+2)*(k+1),2*j*(48+2*ell-2*j-2))
        z+=coeff*(w*m)**j*a**k
    return z


def affine_counts(w):
    m=(6,8,10)
    count=ordinary_counts()
    assert count==[1,0,0,52416000,39007332000,6609020221440]
    n=[min(w//2+(a-6)//2,isqrt(w*a))+1 for a in m]
    ix=(0,n[0],n[0]+n[1])
    sz=sum(n)
    T=[]
    for i,mi in enumerate(m):
        value=Q(count[mi//2])
        for j in range(6):
            row=[0]*(sz+1)
            for a in range(n[i]):
                row[ix[i]+a]=a**(2*j)
            row[-1]=value
            T.append(row)
            value*=Q(w*mi*(2*j+1),48+2*j)
    for ell in range(12,35,2):
        rows=[[zonal(a,w,mi,ell) for a in range(n[i])]
              for i,mi in enumerate(m)]
        for r in Fourier_relations(ell-12):
            T.append([r[i]*x for i,row in enumerate(rows) for x in row]+[0])
    for i in range(1,3):
        for j in range(i):
            a=w//2+(m[i]-m[j])//2
            b=w-a
            if a<n[i] and 0<b<n[j]:
                row=[0]*(sz+1)
                row[ix[i]+a],row[ix[j]+b]=1,-1
                T.append(row)
    piv,R=rref(T,sz)
    free=[j for j in range(sz) if j not in piv]
    B=[None]*sz
    for i,p in enumerate(piv):
        B[p]=tuple(-R[i][j] for j in free)+(R[i][-1],)
    for i,j in enumerate(free):
        B[j]=tuple(Q(k==i) for k in range(len(free)+1))
    return n,ix,B


def small_det(A):
    n=len(A)
    if n==0:
        return 1
    if n==1:
        return A[0][0]
    if n==2:
        return A[0][0]*A[1][1]-A[0][1]*A[1][0]
    return sum((-1)**j*a*small_det([x[:j]+x[j+1:] for x in A[1:]])
               for j,a in enumerate(A[0]))


def lower_counts(B):
    d=len(B[0])-1
    out=[None]*len(B)
    if d==0:
        return [b[-1] for b in B]
    ints=[]
    for row in B:
        l=lcm(*(x.denominator for x in row))
        z=tuple(int(l*x) for x in row)
        common=gcd(*z)
        if common:
            ints.append(tuple(x//common for x in z))
    normals=tuple(r for r in sorted(set(ints)) if any(r[:-1]))
    checked=0
    for rows in combinations(normals,d):
        den=small_det([r[:-1] for r in rows])
        if den==0:
            continue
        nums=[(-1 if (j-d)&1 else 1) *
              small_det([r[:j]+r[j+1:] for r in rows]) for j in range(d)]
        if den<0:
            den=-den
            nums=[-a for a in nums]
        if any(sum(a*x for a,x in zip(r,nums))+r[-1]*den<0 for r in ints):
            continue
        checked+=1
        for i,r in enumerate(B):
            v=(sum(a*x for a,x in zip(r,nums))+r[-1]*den)/den
            assert v>=0
            if out[i] is None or v<out[i]:
                out[i]=v
    assert checked>0 and all(v is not None for v in out)
    return out


EXPECTED={
  12:[72,71,68,63,56,71,60],
  14:[84,83,80,75,68,59,76,63],
  16:[96,95,92,87,80,71,92,79,64],
  18:[108,107,104,99,92,83,72,95,80,99],
  20:[120,119,116,111,104,95,84,111,96,119,100]}


def check(w):
    n,ix,B=affine_counts(w)
    assert len(B[0])-1=={12:1,14:1,16:2,18:3,20:4}[w]
    v=lower_counts(B)
    assert v[0]>0
    if w==18:
        assert v[6]==2256 and v[ix[2]+9]==3688740
    if w==20:
        assert v[ix[2]+9]==Q(242003200,13) and v[ix[2]+10]==299928
    D=[6*w]*(w//2+1)
    for a in range(1,w//2+1):
        for i,m in enumerate((6,8,10)):
            if a<n[i] and v[ix[i]+a]>0:
                D[a]=min(D[a],m*w-a*a)
    if w==16:
        D[8]=min(D[8],4*w)  # the exact norm-4 neighbor lemma
    assert D==EXPECTED[w]
    def d(i,j):
        a=abs(i-j)
        return D[min(a,w-a)]
    best=Q(0)
    for t in sorted(set(D)):
        def cap(A,size):
            nonlocal top
            top=max(top,size)
            while len(A)+size>top:
                vertex,*A=A
                cap([y for y in A if d(vertex,y)>=t],size+1)
        top=1
        cap([i for i in range(1,w) if d(0,i)>=t],1)
        value=Q(top**2,w)*Q(t,4*w)**47
        assert value<Q(47**47,2**236)
        best=max(best,value)
    assert best==Q(4,w)*Q(6*w-1,4*w)**47
    print(w,'PASS, clearance numerators',D,'density upper',float(best)**0.5)


if __name__=='__main__':
    for w in list(map(int,sys.argv[1:])) or EXPECTED:
        check(w)
